/**
 * @ file
 * @ brief NPDE homework about Sobolev evolution problems
 * @ author Ralf Hiptmair
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace SobolevEVP {
/**
 * @brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin stiffness and mass matrices
 *
 * This function first annihilates all selected rows and columns of a matrix in
 * triplet format. Then the corresponding diagonal entries are set to 1, thus
 * preserving the invertibilty of the matrix
 *
 * @param selectvals The predicate identifying selecting the rows and columns to
 * be set to zero
 * @param A matrix in LehrFEM++ internal triplet format. Will be modified!
 */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsColumns(SELECTOR&& selectvals,
                           lf::assemble::COOMatrix<SCALAR>& A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  // Set the selected rows and columns to zero
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  // Set the diagonal entries of zeroed out rows and columns to 1
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/**
 * @brief Assemble Galerkin matrix for scalar linear pure diffusion problem,
 * piecewise linear finite-element space, and homogeneous Dirichlet boundary
 * conditions.
 *
 * @tparam MESHFUNCTION a LehrFEM++ MeshFunction compatible type.
 * @param fe_space_p A description of the used Lagrangian finite element space
 * @param mf_coeff A MeshFunction object providing the diffusion coefficient.
 */
template <typename SCALAR, typename MESHFUNCTION>
Eigen::SparseMatrix<SCALAR> getFEMatrixDirichlet(
    std::shared_ptr<const lf::fe::ScalarFESpace<SCALAR>> fe_space_p,
    const MESHFUNCTION& mf_coeff) {
  // Step I: Building of full Galerkin matrix
  // Obtain local-to-global index mapper
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Provider object for element matrices
  lf::fe::DiffusionElementMatrixProvider<double, MESHFUNCTION> elmat_builder(
      fe_space_p, mf_coeff);
  // Object for sparse matrix to be filled by cell-oriented assembly
  const int N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Invoke cell-oriented assembly of the finite-element Galerkin matrix
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Step II: Take into account Dirichlet boundary conditions
  // The underlying finite-element mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{fe_space_p->Mesh()};
  // Obtain predicate selecting edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_ed_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Boundary flags for actual degrees of freedom
  std::vector<bool> bd_dof_flags(N_dofs, false);
  // Visit all edges of the mesh, retrieve associated dofs and mark them as
  // lying on the boundary
  using gdof_idx_t = lf::assemble::gdof_idx_t;
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    if (bd_ed_flags(*edge)) {
      // Fetch all dof indices associated with the current edge
      std::span<const gdof_idx_t> ed_dof_idx{dofh.GlobalDofIndices(*edge)};
      for (const gdof_idx_t dof_idx : ed_dof_idx) {
        LF_ASSERT_MSG(dof_idx < dofh.NumDofs(),
                      "Dof idx exceeds vector length!");
        bd_dof_flags[dof_idx] = true;
      }
    }
  }
  // Finally set to zero all non-diagonal entries associated with dofs on the
  // boundary
  dropMatrixRowsColumns(
      [&bd_dof_flags](gdof_idx_t idx) -> bool { return bd_dof_flags[idx]; }, A);
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  return A_crs;
}

/**
 * @brief Explicit RK-SSM for MOL for Sobolev evolution problem
 *
 * @tparam MESHFUNCTION_BETA mesh function type for passing coefficient beta
 * @tparam MESHFUNCTION_ALPHA mesh function type for passing coefficient alpha
 * @param fe_space_p pointer to underlying Lagrangian finite element space
 * @param beta coefficient beta wrapper into a MeshFunction object
 * @param alpha coefficient beta wrapper into a MeshFunction object
 * @param mu0 coefficient vector for initial data
 * @param T final time
 * @param A Butcher matrix for explicit RK-SSM, only strict lower triangle is
 * used.
 * @param b weight vector for RK-SSM
 * @M number of timesteps
 */

/* SAM_LISTING_BEGIN_7 */
template <typename MESHFUNCTION_BETA, typename MESHFUNCTION_ALPHA>
Eigen::VectorXd solveRKSobEvl(
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fe_space_p,
    const MESHFUNCTION_BETA& beta, const MESHFUNCTION_ALPHA& alpha,
    const Eigen::VectorXd& mu0, double T, const Eigen::MatrixXd& RK_Mat,
    const Eigen::VectorXd& b, unsigned int M) {
  LF_ASSERT_MSG(mu0.size() == (fe_space_p->LocGlobMap()).NumDofs(),
                "Wrong length of coefficient vector");
  // Build Galerkin matrices taking into account homogeneous Dirichlet boundary
  // conditions
  Eigen::SparseMatrix<double> B =
      getFEMatrixDirichlet<double, MESHFUNCTION_BETA>(fe_space_p, beta);
  Eigen::SparseMatrix<double> A =
      getFEMatrixDirichlet<double, MESHFUNCTION_ALPHA>(fe_space_p, alpha);
  Eigen::VectorXd muj(mu0);  // current state vector
  // LU-decompose B outside timestepping loop, which is important for efficient
  // implementation
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(B);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");

  int s = RK_Mat.cols();  // Number of stages
  LF_ASSERT_MSG(s == RK_Mat.rows(), "Butcher matrix must be square!");
  LF_ASSERT_MSG(s == b.size(), "s weights required!");
  // timestep size
  double tau = T / M;
  // Main timestepping loop
  Eigen::VectorXd mu_next(mu0.size());            // next state vector
  Eigen::VectorXd tmp(mu0.size());                // Summation vector
  std::vector<Eigen::VectorXd> incs(s, mu_next);  // RK increments
  // The r.h.s. vector field for MOL ODE in standard form
  // Invoking solver.solve() amounts to the application of the
  // inverse of the matrix B.
  auto Vf = [&solver, &A](const Eigen::VectorXd& y) -> Eigen::VectorXd {
    const Eigen::VectorXd fy = -solver.solve(A * y);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
    return fy;
  };
  for (int k = 0; k < M; ++k) {
    // First increment vector
    incs[0] = Vf(muj);
    // Compute increments successively, which is possible for an explicit RK-SSM
    // Simultaneous assemble next state vector
    mu_next = muj + (tau * b[0]) * incs[0];
    for (int i = 1; i < s; ++i) {
      // Weighted sum of already computed increments
      tmp.setZero();
      for (int j = 0; j < i; ++j) {
        tmp += RK_Mat(i, j) * incs[j];
      }
      incs[i] = Vf(muj + tau * tmp);
      // Update state by weighed sum of increments
      mu_next += (tau * b[i]) * incs[i];
    }
    // Efficient way to set current state to new state
    std::swap(muj, mu_next);
  }
  return muj;
}
/* SAM_LISTING_END_7 */

}  // namespace SobolevEVP
