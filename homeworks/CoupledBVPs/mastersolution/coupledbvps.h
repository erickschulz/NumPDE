/**
 * @file coupledbvps.h
 * @brief NPDE homework CoupledBVPs code
 * @author W. Tonnon
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/fe/fe.h>
#include <lf/io/gmsh_reader.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <memory>

namespace CoupledBVPs {

// Solve Laplace equation with coefficient alpha, source term f, and dirichlet
// boundary term g
/* SAM_LISTING_BEGIN_7 */
template <typename MESHFN_ALPHA, typename MESHFN_F, typename G_FUNCTOR>
Eigen::VectorXd solveDirichletBVP(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    MESHFN_ALPHA alpha, MESHFN_F F, G_FUNCTOR&& g) {
  // Obtain the mesh
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Dimension of finite element space`
  const lf::base::size_type N_dofs(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  // Define local computations for the LHS through a matrix provider
  lf::fe::DiffusionElementMatrixProvider elmat_builder(fes_p, alpha);
  // Construct the full matrix through local computations
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Define local computations for the RHS through a vector provider
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(F)>
      elvec_builder(fes_p, F);

  // Define the RHS vector
  Eigen::VectorXd phi(N_dofs);
  // We initialize to zero, because local assembly adds to the vector
  phi.setConstant(0);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Flag all nodes on the boundary (and only those)
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Set up predicate: Run through all global shape functions and check whether
  // they are associated with an entity on the boundary, store Dirichlet data.
  std::vector<std::pair<bool, double>> ess_dof_select{};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
    const lf::mesh::Entity& dof_node{dofh.Entity(dofnum)};
    const Eigen::Vector2d node_pos{
        lf::geometry::Corners(*dof_node.Geometry()).col(0)};
    const double g_val = g(node_pos);
    if (bd_flags(dof_node)) {
      // Dof associated with a entity on the boundary: "essential dof"
      // The value of the dof should be set to the value of the function
      // u at the location of the node.
      ess_dof_select.emplace_back(true, g_val);
    } else {
      // Interior node, also store value of solution for comparison purposes
      ess_dof_select.emplace_back(false, g_val);
    }
  }
  // modify linear system of equations
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&ess_dof_select](int dof_idx) -> std::pair<bool, double> {
        return ess_dof_select[dof_idx];
      },
      A, phi);

  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  // Return the solution
  return sol_vec;
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
template <typename F_FUNCTOR>
std::pair<Eigen::VectorXd, Eigen::VectorXd> solveModulatedHeatFlow(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    F_FUNCTOR&& f) {
  // We need the following lambda functions to call solveDirichletBVP
  auto f2 = [&f](Eigen::Vector2d x) -> double { return f(x) * f(x); };
  auto g = [](Eigen::Vector2d /*x*/) -> double { return 1.; };
  auto zero = [](Eigen::Vector2d /*x*/) -> double { return 0.; };
  // We need the following MeshFunctions to call solveDirichletBVP
  lf::mesh::utils::MeshFunctionGlobal mf_f2{f2};
  lf::mesh::utils::MeshFunctionConstant mf_one(1.);

  Eigen::VectorXd w_gf;  // Basis expansion coefficients for $w_h$
  Eigen::VectorXd v_gf;  // Basis expansion coefficients for $v_h$
  // We first solve for w using the provided solver for diffusion problems
  // (solveDirichletBVP)
  w_gf = solveDirichletBVP(fes_p, mf_one, mf_f2, g);
  // We cast w into a meshfunction such that it can be used in solveDirichletBVP
  // as a coefficient and source
  lf::fe::MeshFunctionFE mf_w(fes_p, w_gf);
  // We solve for v using solveDirichletBVP using the previously computed source
  // and coefficient w
  v_gf = solveDirichletBVP(fes_p, mf_w, mf_w, zero);
  // We return the solution
  return std::make_pair(w_gf, v_gf);
}
/* SAM_LISTING_END_8 */

}  // namespace CoupledBVPs
