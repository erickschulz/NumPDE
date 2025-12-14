/**
 * @file taylorhoodnonmonolithic.h
 * @brief NPDE homework TaylorHoodNonMonolithic code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef TaylorHoodNonMonolithic_H_
#define TaylorHoodNonMonolithic_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/assembly_types.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/fe_space_lagrange_o2.h>
#include <lf/uscalfe/precomputed_scalar_reference_finite_element.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <memory>

namespace TaylorHoodNonMonolithic {

/** @brief
 * ENTITY_MATRIX_PROVIDER for bilinear form b
 *
 */
enum DirFlag : Eigen::Index { X_Dir = 0, Y_Dir = 1 };
class THBElementMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, 3, 6>;
  // Main constructor
  THBElementMatrixProvider(const THBElementMatrixProvider&) = delete;
  THBElementMatrixProvider(THBElementMatrixProvider&&) noexcept = default;
  THBElementMatrixProvider& operator=(const THBElementMatrixProvider&) = delete;
  THBElementMatrixProvider& operator=(THBElementMatrixProvider&&) = delete;
  THBElementMatrixProvider(DirFlag dirflag) : dirflag_(dirflag) {}
  virtual ~THBElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);

 private:
  DirFlag dirflag_;  // Chooses partial derivative in bilinear form
};

/** @brief Assembly of (extended) saddle point linear system
 *
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FFUNCTOR>
std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>,
           Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::VectorXd>
buildStokesLSE(std::shared_ptr<const lf::mesh::Mesh> mesh_p, FFUNCTOR&& force) {
  LF_ASSERT_MSG(mesh_p != nullptr, "Mesh must be supplied!");
  // I. Assemble full Galerkin matrix A in triplet format
  // Initialize quadratic Lagrange FE space
  auto fes_velo_comp =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Unit "coefficient function", $\mu=1$!
  auto one = [](const Eigen::Vector2d& /*x*/) -> double { return 1.0; };
  // Create unit mesh function
  const lf::mesh::utils::MeshFunctionGlobal mf_one{one};
  // Element matrix provider to $-\Delta$ bilinear form
  lf::fe::DiffusionElementMatrixProvider Laplace_emp(fes_velo_comp, mf_one);
  // local $\to$ global index mapping
  const lf::assemble::DofHandler& dofh_velo_comp{fes_velo_comp->LocGlobMap()};
  const lf::assemble::size_type N_velo_comp = dofh_velo_comp.NumDofs();
  // Full matrix A in triplet format
  lf::assemble::COOMatrix<double> A_velo_comp(N_velo_comp, N_velo_comp);
  lf::assemble::AssembleMatrixLocally(0, dofh_velo_comp, dofh_velo_comp,
                                      Laplace_emp, A_velo_comp);
  /* SAM_LISTING_END_1 */
  /* SAM_LISTING_BEGIN_2 */
  // II. Assemble right hand side vectors $\vec{\varphibf}_x$,
  // $\vec{\varphibf}_y$
  // Set up mesh functions for the two components of the force field
  auto f_x = [&force](Eigen::Vector2d x) -> double { return force(x)[0]; };
  auto f_y = [&force](Eigen::Vector2d x) -> double { return force(x)[1]; };
  lf::mesh::utils::MeshFunctionGlobal mf_f_x(f_x);
  lf::mesh::utils::MeshFunctionGlobal mf_f_y(f_y);
  // Assembly of full right-hand side vectors
  lf::uscalfe::ScalarLoadElementVectorProvider phi_x_builder(fes_velo_comp,
                                                             mf_f_x);
  lf::uscalfe::ScalarLoadElementVectorProvider phi_y_builder(fes_velo_comp,
                                                             mf_f_y);
  Eigen::VectorXd phi_x(N_velo_comp);
  phi_x.setZero();
  Eigen::VectorXd phi_y(N_velo_comp);
  phi_y.setZero();
  // Invoke assembly on cells (codim == 0)
  lf::assemble::AssembleVectorLocally(0, dofh_velo_comp, phi_x_builder, phi_x);
  lf::assemble::AssembleVectorLocally(0, dofh_velo_comp, phi_y_builder, phi_y);

  // III. Take into account {\bf homogeneous} essential boundary conditions
  // Flag \cor{any} entity located on the boundary
  const auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
  // Flag d.o.f.s on the boundary
  std::vector<bool> dof_bd_flags(N_velo_comp, false);
  // The trouble is that FixFlaggedSolutionComponents() does not allow to pass
  // multiple right-hand side vectors. Therefore we have to set the components
  // of phi\_y corrresponding to d.o.f.s on the boundary to zero "manually".
  for (lf::assemble::glb_idx_t dof_idx = 0; dof_idx < N_velo_comp; ++dof_idx) {
    const bool on_bd = bd_flags(dofh_velo_comp.Entity(dof_idx));
    dof_bd_flags[dof_idx] = on_bd;
    if (on_bd) {
      phi_y[dof_idx] = 0.0;
    }
  }
  // Modify linear system of equations and r.h.s. vector phi\_x
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&dof_bd_flags](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return {dof_bd_flags[dof_idx], 0.0}; },
      A_velo_comp, phi_x);
  // IV. Aseemble the matrices $\VB_x$ and $\VB_y$
  // Initialize lowest-order Lagrange FE space, p=1
  auto fes_pressure =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // local $\to$ global index mapping and number of d.o.f.s
  const lf::assemble::DofHandler& dofh_pressure{fes_pressure->LocGlobMap()};
  const lf::assemble::size_type n_pressure = dofh_pressure.NumDofs();
  // Helper objects for computation of element matrices
  THBElementMatrixProvider B_emp_x(X_Dir);
  THBElementMatrixProvider B_emp_y(Y_Dir);
  // Matrices in triplet format
  lf::assemble::COOMatrix<double> B_x(n_pressure, N_velo_comp);
  lf::assemble::COOMatrix<double> B_y(n_pressure, N_velo_comp);
  lf::assemble::AssembleMatrixLocally(0, dofh_velo_comp, dofh_pressure, B_emp_x,
                                      B_x);
  lf::assemble::AssembleMatrixLocally(0, dofh_velo_comp, dofh_pressure, B_emp_y,
                                      B_y);

  // V. Set columns of B\_x and B\_y corresponding to d.o.f. on the boundary to
  // zero.
  auto zero_col_sel = [&dof_bd_flags](
                          lf::assemble::gdof_idx_t /*idx_row*/,
                          lf::assemble::gdof_idx_t idx_col) -> bool {
    return dof_bd_flags[idx_col];
  };
  B_x.setZero(zero_col_sel);
  B_y.setZero(zero_col_sel);

  return {A_velo_comp.makeSparse(), B_x.makeSparse(), B_y.makeSparse(), phi_x,
          phi_y};
}
#endif
/* SAM_LISTING_END_2 */

/** @brief CG Uzawa algorithm for solving saddle-point linear system of
 * equations
 *
 */

/* SAM_LISTING_BEGIN_5 */
template <typename RECORDER =
              std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&,
                                 const Eigen::VectorXd&)>>
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> CGUzawa(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B_x,
    const Eigen::SparseMatrix<double>& B_y, const Eigen::VectorXd& phi_x,
    const Eigen::VectorXd& phi_y, double rtol = 1E-6, double atol = 1E-8,
    unsigned int itmax = 100,
    RECORDER&& rec = [](const Eigen::VectorXd& mu_x,
                        const Eigen::VectorXd& mu_y,
                        const Eigen::VectorXd& pi) -> void {}) {
  const Eigen::Index N = A.cols();
  const Eigen::Index M = B_x.rows();
  // Check consistent sizes of matrices
  LF_VERIFY_MSG(A.cols() == A.rows(), "A must be square");
  LF_VERIFY_MSG(B_x.cols() == N, "B_x.cols() == A.cols() required");
  LF_VERIFY_MSG(B_y.cols() == N, "B_y.cols() == A.cols() required");
  LF_VERIFY_MSG(B_y.rows() == M, "B_x and B_y must have the same size");
  // (Sparse) LU-decomposition of A: done once for the sake of efficiency!
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
  solver_A.compute(A);
  LF_VERIFY_MSG(solver_A.info() == Eigen::Success, "LU decomposition failed");
  // Initialize the remaining vectors
  Eigen::VectorXd pi(M);
  pi.setZero();
  Eigen::VectorXd mu_x = solver_A.solve(phi_x);
  LF_VERIFY_MSG(solver_A.info() == Eigen::Success, "Solving x LSE failed");
  Eigen::VectorXd mu_y = solver_A.solve(phi_y);
  LF_VERIFY_MSG(solver_A.info() == Eigen::Success, "Solving y LSE failed");
  Eigen::VectorXd delta = B_x * mu_x + B_y * mu_y;
  Eigen::VectorXd rho = -delta;
  const double rn0 = rho.norm();
  double rn = rn0;
  Eigen::VectorXd omega_x(N);
  Eigen::VectorXd omega_y(N);
  Eigen::VectorXd eta_x(N);
  Eigen::VectorXd eta_y(N);
  // Main iteration loop
  unsigned int step = 0;
  do {
    if (rn >= atol) {
      //  $\cob{\vec{\omegabf} := \cod{\VB^{\top}}\vec{\deltabf}}$;
      omega_x = B_x.transpose() * delta;
      omega_y = B_y.transpose() * delta;
      //  $\cob{\vec{\etabf} := \cop{\VA^{-1}}\vec{\omegabf}}$;
      eta_x = solver_A.solve(omega_x);
      eta_y = solver_A.solve(omega_y);
      // $\cob{\alpha}$ = (rn * rn) / $\cob{\vec{\omegabf}^{\top}\vec{\etabf}}$;
      const double alpha =
          (rn * rn) / (omega_x.dot(eta_x) + omega_y.dot(eta_y));
      // $\cob{\vec{\pibf}:= \vec{\pibf} + \alpha \vec{\deltabf}}$;
      pi += alpha * delta;
      // $\cob{\vec{\mubf} := \vec{\mubf} - \alpha \vec{\etabf}}$;
      mu_x -= alpha * eta_x;
      mu_y -= alpha * eta_y;
      // $\cob{\vec{\rhobf} := \vec{\gammabf} - \cod{\VB}\vec{\mubf}}$;
      rho = -(B_x * mu_x + B_y * mu_y);
      const double rn1 = rho.norm();
      const double beta = (rn1 * rn1) / (rn * rn);
      rn = rn1;
      // $\cob{\vec{\deltabf} := -\vec{\rhobf} + \beta \vec{\deltabf}}$;
      delta = -rho + beta * delta;
    }
    rec(mu_x, mu_y, pi);
    step++;
  } while ((rn >= rn0 * rtol) and (rn >= atol) and (step < itmax));
  return {mu_x, mu_y, pi};
}
/* SAM_LISTING_END_5 */

/** @brief Study convergence of Uzawa iteration
 *
 */
Eigen::Matrix<double, Eigen::Dynamic, 6> monitorUzawaConvergence(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, bool print = false);

}  // namespace TaylorHoodNonMonolithic
