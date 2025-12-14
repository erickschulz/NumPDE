/**
 * @file stokespipeflow.h
 * @brief NPDE homework StokesPipeFlow code
 * @ author Wouter Tonnon
 * @ date May 2024
 * @ copyright Developed at SAM, ETH Zurich
 */

#ifndef StokesPipeFlow_H_
#define StokesPipeFlow_H_

// Include almost all parts of LehrFEM++; some my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembly_types.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/lf_assert.h>
#include <lf/base/ref_el.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/special_entity_sets.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cstddef>

namespace StokesPipeFlow {
/**
 * @brief Element matrix provider for Taylor-Hood Stokes FEM
 */
/* SAM_LISTING_BEGIN_1 */
class TaylorHoodElementMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, 15, 15>;
  TaylorHoodElementMatrixProvider(const TaylorHoodElementMatrixProvider&) =
      delete;
  TaylorHoodElementMatrixProvider(TaylorHoodElementMatrixProvider&&) noexcept =
      default;
  TaylorHoodElementMatrixProvider& operator=(
      const TaylorHoodElementMatrixProvider&) = delete;
  TaylorHoodElementMatrixProvider& operator=(
      TaylorHoodElementMatrixProvider&&) = delete;
  TaylorHoodElementMatrixProvider() = default;
  virtual ~TaylorHoodElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);

 private:
  ElemMat AK_;
};
/* SAM_LISTING_END_1 */

/**
 * @brief Assembly of full Galerkin matrix in triplet format
 *
 * @param dofh DofHandler object for all FE spacesw
 */
lf::assemble::COOMatrix<double> buildTaylorHoodGalerkinMatrix(
    const lf::assemble::DofHandler& dofh);

/**
 * @brief Taylor-Hood FE solultion of pipe flow problem
 *
 * @tparam functor type taking a 2-vector and returning a 2-vector
 * @param dofh DofHandler object for all FE spaces
 * @param g functor providing Dirchlet boundary data
 */
/* SAM_LISTING_BEGIN_2 */
template <typename gFunctor>
Eigen::VectorXd solvePipeFlow(const lf::assemble::DofHandler& dofh,
                              gFunctor&& g, bool print = true) {
  // Number of d.o.f. in FE spaces
  size_t n = dofh.NumDofs();
  if (print)
    std::cout << "Computing: solvePipeFlow() with N = " << n
              << ", assembling .. " << std::flush;
  // Obtain full Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A{buildTaylorHoodGalerkinMatrix(dofh)};
  if (print) std::cout << "done. Imposing BDC ... " << std::flush;
  LF_VERIFY_MSG(A.cols() == A.rows(), "Matrix A must be square");
  // Auxiliary right-hnad side vector
  Eigen::VectorXd phi(A.cols());
  phi.setZero();
  // Impose Dirichlet boundary conditions
#if SOLUTION
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  // Flag \cor{any} entity located on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
  // Flag vector for d.o.f. on the boundary
  std::vector<std::pair<bool, double>> ess_dof_select(n + 1, {false, 0.0});
  // Visit nodes on the boundary
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    if (bd_flags(*node)) {
      // Indices of global shape functions sitting at node
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*node)};
      LF_ASSERT_MSG(dof_idx.size() == 3, "Node must carry 3 dofs!");
      // Position of node
      const Eigen::Vector2d pos{Corners(*(node->Geometry())).col(0)};
      // Dirichlet data
      const Eigen::Vector2d g_val{g(pos)};
      // x-component of the velocity
      ess_dof_select[dof_idx[0]] = {true, g_val[0]};
      // y-component of the velocity
      ess_dof_select[dof_idx[1]] = {true, g_val[1]};
    }
  }
  // Visit edges on the boundasry
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    if (bd_flags(*edge)) {
      // Indices of global shape functions associated with the edge
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*edge)};
      LF_ASSERT_MSG(dof_idx.size() == 2, "Edge must carry 2 dofs!");
      // Midpoint of edge
      const Eigen::MatrixXd endpoints{Corners(*(edge->Geometry()))};
      const Eigen::Vector2d pos{0.5 * (endpoints.col(0) + endpoints.col(1))};
      // Dirichlet data
      const Eigen::Vector2d g_val{g(pos)};
      // x-component of the velocity
      ess_dof_select[dof_idx[0]] = {true, g_val[0]};
      // y-component of the velocity
      ess_dof_select[dof_idx[1]] = {true, g_val[1]};
    }
  }
  // modify linear system of equations
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      A, phi);
#else
  /* **********************************************************************
     Your code here
     ********************************************************************** */
#endif
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  if (print) std::cout << "done. Solving ..... " << std::flush;
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  // Examine return status of solver in case the matrix is singular
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  const Eigen::VectorXd dofvec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  if (print) std::cout << "done. |dof vector| = " << dofvec.norm() << std::endl;
  // This is the coefficient vector for the FE solution; Dirichlet
  // boundary conditions are included
  return dofvec;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Compute dissipated power for a Talor-Hood FE solution
 *
 * @param dofh DofHandler object for monolithic Taylor-Hood FEM
 * @param mu_vec basis expansion coefficient vector
 */
double compDissPowVolume(const lf::assemble::DofHandler& dofh,
                         const Eigen::VectorXd& mu_vec);

/**
 * @brief Computes dissipated power by boundary-based formula
 *
 * @param dofh DofHandler object for monolithic Taylor-Hood FEM
 * @param mu_vec basis expansion coefficient vector
 *
 * This implementation is valid only in the special geometric setting of the
 * pipe flow model for the homework project StokesPipeFLow.
 */
double compDissPowBd(const lf::assemble::DofHandler& dofh,
                     const Eigen::VectorXd& muvec, bool print = false);

/**
 * @brief Convergence test for Tyalor-Hood FEM
 */
void testCvgTaylorHood(unsigned int refsteps = 5);

enum PowerFlag { NOCMOP, VOLUME, BOUNDARY };

/**
 * @brief Taylor-Hood FEM for pipe flow: visualization of flow field and
 *computation of dissipated power by different formulas.
 *
 * @param powerflag: If NOCMOP, do not compute dissipated power, if VOLUME, use
 * integral of squared norm of curl v, if BOUNDARY use pressure-based integral
 * over boundary.
 * @param producevtk: If true write solution to .vtk file for visualization with
 * Paraview
 * @param meshfile name of the gmsh mesh file to read triangulation from
 * @param outfile base name of .vtk output files
 *
 */
double allPipeFlow(PowerFlag powerflag, bool producevtk, const char* meshfile,
                   const char* outfile = nullptr);

/**
 * @brief Visualization of FEM solution for pipe flow setting
 *
 * @param meshfile name of the gmsh mesh file to read triangulation from
 * @param outfile base name of .vtk output files
 */
void visualizeTHPipeFlow(const char* meshfile = "pipe.msh",
                         const char* outfile = "pipeflow");

/**
 * @brief Compute dissipated power based on Taylor-Hood FEM simulation:
 * volume-integration based formuls
 */
double computeDissipatedPower(const char* meshfile = "pipe.msh");
/**
 * @brief Compute dissipated power based on Taylor-Hood FEM simulation:
 * volume-integration based formuls
 */
double computeDissipatedPowerBd(const char* meshfile = "pipe.msh");

}  // namespace StokesPipeFlow

#endif
