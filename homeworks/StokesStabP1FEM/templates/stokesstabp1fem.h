/**
 * @file stokesstabp1fem.h
 * @brief NPDE homework StokesStabP1FEM code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef StokesStabP1FEM_H_
#define StokesStabP1FEM_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembler.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace StokesStabP1FEM {
/**
 * @brief Element matrix provider for Taylor-Hood Stokes FEM
 */
/* SAM_LISTING_BEGIN_1 */
class P1StabFEMElementMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, 9, 9>;
  P1StabFEMElementMatrixProvider(const P1StabFEMElementMatrixProvider&) =
      delete;
  P1StabFEMElementMatrixProvider(P1StabFEMElementMatrixProvider&&) noexcept =
      default;
  P1StabFEMElementMatrixProvider& operator=(
      const P1StabFEMElementMatrixProvider&) = delete;
  P1StabFEMElementMatrixProvider& operator=(P1StabFEMElementMatrixProvider&&) =
      delete;
  P1StabFEMElementMatrixProvider() = default;
  virtual ~P1StabFEMElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);

 private:
  ElemMat MK_;
};
/* SAM_LISTING_END_1 */

/**
 * @brief ENTITY_VECTOR_PROVIDER class for r.h.s. forcing term
 * for P1-FEM for Stokes
 *
 */
/* SAM_LISTING_BEGIN_7 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
class P1StabFEMElementVectorProvider {
 public:
  using ElemVec = Eigen::Matrix<double, 9, 1>;
  P1StabFEMElementVectorProvider(const P1StabFEMElementVectorProvider&) =
      delete;
  P1StabFEMElementVectorProvider(P1StabFEMElementVectorProvider&&) noexcept =
      default;
  P1StabFEMElementVectorProvider& operator=(
      const P1StabFEMElementVectorProvider&) = delete;
  P1StabFEMElementVectorProvider& operator=(P1StabFEMElementVectorProvider&&) =
      delete;
  virtual ~P1StabFEMElementVectorProvider() = default;

  P1StabFEMElementVectorProvider(MESH_FUNCTION f) : f_(f) {}
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemVec Eval(const lf::mesh::Entity& cell);

 private:
  ElemVec phiK_;
  MESH_FUNCTION f_;
};
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
typename P1StabFEMElementVectorProvider<MESH_FUNCTION>::ElemVec
P1StabFEMElementVectorProvider<MESH_FUNCTION>::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Obtain values of forcing vector field in midpoints of edges
  // Matrix collecting the \textbf{reference coordinates} of the midpoints of
  // the edges of the triangle
  // clang-format off
    const Eigen::MatrixXd mp{(Eigen::MatrixXd(2, 3) <<
			      0.5, 0.5, 0.0,
			      0.0, 0.5, 0.5).finished()};
  // clang-format on
  std::vector<Eigen::Vector2d> f_vals{f_(cell, mp)};
  LF_ASSERT_MSG(f_vals.size() == 3, "Expect 3 f values");
  phiK_.setZero();
  phiK_.segment<2>(0) = area / 6.0 * (f_vals[0] + f_vals[2]);
  phiK_.segment<2>(3) = area / 6.0 * (f_vals[0] + f_vals[1]);
  phiK_.segment<2>(6) = area / 6.0 * (f_vals[1] + f_vals[2]);
  return phiK_;
}
/* SAM_LISTING_END_8 */

/**
 * @brief Assembly of full Galerkin matrix in triplet format
 *
 * @param dofh DofHandler object for all FE spaces (= monolithic FE space)
 */
lf::assemble::COOMatrix<double> buildP1StabFEMGalerkinMatrix(
    const lf::assemble::DofHandler& dofh);

/**
 * @brief Taylor-Hood FE solultion of pipe flow problem
 *
 * @tparam functor type taking a 2-vector and returning a 2-vector
 * @param dofh DofHandler object for all FE spaces
 * @param g functor providing Dirchlet boundary data
 */
/* SAM_LISTING_BEGIN_2 */
template <lf::mesh::utils::MeshFunction FUNCTOR_F>
Eigen::VectorXd solveP1StabFEMStokesBVP(const lf::assemble::DofHandler& dofh,
                                        FUNCTOR_F f, bool print = true) {
  // Number of d.o.f. in FE spaces
  size_t n = dofh.NumDofs();
  if (print)
    std::cout << "Computing: solveP1StabFEMStokesBVP( with N = " << n
              << ", assembling .. " << std::flush;
  // Step I: Building full linear system of equatins
  // Obtain full Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A{buildP1StabFEMGalerkinMatrix(dofh)};
  if (print) std::cout << "done. Imposing BDC ... " << std::flush;
  LF_VERIFY_MSG(A.cols() == A.rows(), "Matrix A must be square");
  // Build right-hnad side vector
  Eigen::VectorXd phi(A.cols());
  phi.setZero();
  P1StabFEMElementVectorProvider f_evp(f);
  lf::assemble::AssembleVectorLocally(0, dofh, f_evp, phi);
  // Step II: Impose zero Dirichlet boundary conditions
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  // Flag nodes ocated on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Flag vector for d.o.f. on the boundary
  std::vector<std::pair<bool, double>> ess_dof_select(n + 1, {false, 0.0});
  // Visit nodes on the boundary
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    if (bd_flags(*node)) {
      // Indices of global shape functions sitting at node
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*node)};
      LF_ASSERT_MSG(dof_idx.size() == 3, "Node must carry 3 dofs!");
      // x-component of the velocity
      ess_dof_select[dof_idx[0]] = {true, 0.0};
      // y-component of the velocity
      ess_dof_select[dof_idx[1]] = {true, 0.0};
    }
  }
  // modify linear system of equations
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      A, phi);
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
 * @brief Convergence test for simple Stokes FEM
 */
void testCvgP1StabFEM(unsigned int refsteps);

}  // namespace StokesStabP1FEM

#endif
