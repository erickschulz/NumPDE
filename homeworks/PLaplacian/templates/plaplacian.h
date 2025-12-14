/**
 * @file plaplacian.h
 * @brief NPDE homework PLaplacian code
 * @author W. Tonnon and R. Hiptmair
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <cmath>
#include <iostream>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/fe/loc_comp_ellbvp.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/mesh_function_traits.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/fe_space_lagrange_o1.h>
#include <lf/uscalfe/fe_space_lagrange_o2.h>
#include <lf/uscalfe/uscalfe.h>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace plaplacian {

/* SAM_LISTING_BEGIN_2 */
template <typename MESHFUNCTION, typename FUNCTION>
class FunctionMFWrapper {
 public:
  static_assert(lf::mesh::utils::MeshFunction<MESHFUNCTION>);
  using mf_vector_t = decltype(std::declval<MESHFUNCTION>()(
      std::declval<const lf::mesh::Entity>(),
      std::declval<const Eigen::MatrixXd>()));
  using mf_result_t = typename mf_vector_t::value_type;
  using F_result_t =
      decltype(std::declval<FUNCTION>()(std::declval<mf_result_t>()));

  explicit FunctionMFWrapper(MESHFUNCTION mf, FUNCTION F)
      : mf_(std::move(mf)), F_(std::move(F)) {}
  FunctionMFWrapper(const FunctionMFWrapper&) = default;
  FunctionMFWrapper(FunctionMFWrapper&&) noexcept = default;
  FunctionMFWrapper& operator=(const FunctionMFWrapper&) = delete;
  FunctionMFWrapper& operator=(FunctionMFWrapper&&) = delete;
  ~FunctionMFWrapper() = default;

  std::vector<F_result_t> operator()(const lf::mesh::Entity& e,
                                     const Eigen::MatrixXd& local) const {
    LF_ASSERT_MSG(e.RefEl().Dimension() == local.rows(),
                  "mismatch between entity dimension and local.rows()");
    const std::vector<mf_result_t> mf_result = mf_(e, local);
    std::vector<F_result_t> result(local.cols());
    for (long i = 0; i < local.cols(); ++i) {
      result[i] = F_(mf_result[i]);
    }
    return result;
  }

 private:
  MESHFUNCTION mf_;
  FUNCTION F_;
};
/* SAM_LISTING_END_2 */

/** @brief One step of fixed-point iteration for Dirchlet BVP for p-Laplacian
 *
 * @param fes_p lowest-order Lagrangian FE space (Also contains information
 * about the mesh and the local-to-global index mapping)
 * @param f right-hanmd side source function
 * @param p power of Laplacian
 * @param mu_vec tent function basis expansion coefficient vector of current
 * iterate. This reference parameter will also be used to return the new
 * iterate.
 * @return No return value. Result is returned via mu_vec argument
 */

/* SAM_LISTING_BEGIN_1 */
template <typename F_FUNCTOR>
void fixedPointNextIt(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    F_FUNCTOR f, double p, Eigen::VectorXd& mu_vec) {
  // Wrap right hand side source function into a Mesh Function
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Dimension of finite element space`
  const std::size_t N_dofs(dofh.NumDofs());
  // We will need a vanishing MeshFunction
  lf::mesh::utils::MeshFunctionConstant mf_zero(0.);
  lf::mesh::utils::MeshFunctionConstant mf_one(1.);
  // We need the gradient of the current iteration in MeshFunction format
  const lf::fe::MeshFunctionGradFE mf_grad_uh_prev(fes_p, mu_vec);

  // We need the diffusion coefficient as MeshFunction
  const FunctionMFWrapper mf_abs_uh_pmin2(mf_grad_uh_prev,
                                          [p](Eigen::Vector2d xi) -> double {
// ************************
// Adjust the following code
  return 0.;
// ************************
                                          });

  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  // ************************
  // Adjust the following code
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_one),
                                                      decltype(mf_zero)>
      elmat_provider(fes_p, mf_one, mf_zero);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);
// ************************

  // Enforce homogeneous Dirichlet boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Elimination of degrees of freedom on the boundary. Also sets the
  // corresponding entries of phi to zero.
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};

        const Eigen::Vector2d node_pos{
            lf::geometry::Corners(*node.Geometry()).col(0)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);

  // Convert matrix A in COO format to A\_crs in Eigen::SparseMatrix format
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve the system and store the solution in mu\_vec
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  mu_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return;
}
/* SAM_LISTING_END_1 */

/** @brief One step of Newton iteration for Dirchlet BVP for p-Laplacian
 *
 * @param fes_p lowest-order Lagrangian FE space (Also contains information
 * about the mesh and the local-to-global index mapping)
 * @param f right-hanmd side source function
 * @param p power of Laplacian
 * @param mu_vec tent function basis expansion coefficient vector of current
 * iterate.
 * @return difference between next and current iterate
 */

/* SAM_LISTING_BEGIN_3 */
template <typename F_FUNCTOR>
Eigen::VectorXd newtonUpdate(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    F_FUNCTOR f, double p, const Eigen::VectorXd& mu_vec) {
  // Initialize the return value
  Eigen::VectorXd upd(mu_vec.size());
  // Wrap right hand side source function into a Mesh Function
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Dimension of finite element space`
  const std::size_t N_dofs(dofh.NumDofs());
  // We will need a vanishing MeshFunction
  lf::mesh::utils::MeshFunctionConstant mf_zero(0.);
  // We need the gradient of the current iteration in MeshFunction format
  const lf::fe::MeshFunctionGradFE mf_grad_uh_prev(fes_p, mu_vec);

  // We need ||grad uh\_prev||\^(p-2)*[(p-2)/||grad uh\_prev||\^2*grad
  // uh\_prev*grad uh\_prev\^T + I] as MeshFunction
  const FunctionMFWrapper mf_coeff_LHS(
      mf_grad_uh_prev, [p](Eigen::Vector2d xi) -> Eigen::MatrixXd {
        return std::pow(xi.norm(), p - 2) *
               ((p - 2) / xi.squaredNorm() * xi * xi.transpose() +
                Eigen::Matrix2d::Identity());
      });
  // We need ||grad uh\_prev||\^(p-2) as MeshFunction
  const FunctionMFWrapper mf_coeff_RHS(
      mf_grad_uh_prev,
      [p](Eigen::Vector2d xi) -> double { return std::pow(xi.norm(), p - 2); });

  // Assemble matrix A, representing (||grad uh\_prev||\^(p-2)*[(p-2)/||grad
  // uh\_prev||\^2*grad uh\_prev*grad uh\_prev\^T + I] grad uh, grad vh) in COO
  // format
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_coeff_LHS), decltype(mf_zero)>
      elmat_provider_A(fes_p, mf_coeff_LHS, mf_zero);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider_A, A);

  // Assemble matrix M, representing (||grad uh\_prev||\^(p-2) grad uh, grad vh)
  // in COO format
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_coeff_RHS), decltype(mf_zero)>
      elmat_provider_M(fes_p, mf_coeff_RHS, mf_zero);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider_M, M);

  // Convert M to Eigen::SparseMatrix format
  Eigen::SparseMatrix<double> M_crs = M.makeSparse();

  // Assemble right-hand side vector
  Eigen::VectorXd phi(N_dofs);
  // phi = 0
  phi.setZero();
  // Assemble volume part of right-hand side vector depending on the source
  // function f.
  // Initialize object taking care of local computations on all cells.
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fes_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  // phi = (f,vh)
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // phi = -(||grad uh\_prev||\^(p-2) grad uh\_prev, grad vh)
  phi -= M_crs * mu_vec;

  // Enforce homogeneous Dirichlet boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Elimination of degrees of freedom on the boundary. Also sets the
  // corresponding entries of phi to zero.
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};

        const Eigen::Vector2d node_pos{
            lf::geometry::Corners(*node.Geometry()).col(0)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);

  // Convert matrix A in COO format to A\_crs in Eigen::SparseMatrix format
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve the system and add the solution to mu\_vec
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  upd = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return upd;
}
/* SAM_LISTING_END_3 */

/** @brief Fixed-point-based solver for Dirchlet BVP for p-Laplacian
 *
 * @param fes_p lowest-order Lagrangian FE space (Also contains information
 * about the mesh and the local-to-global index mapping)
 * @param f right-hanmd side source function
 * @param p power of Laplacian
 * @param rtol relative error tolerance
 * @param atol absolute error tolerance
 * @param itmax max fixed-point iterations
 * @param rec recorder to collect data
 * @return tent function basis expansion of computed solution
 */
template <typename FUNCTOR_F, typename RECORDER = std::function<
                                  void(const Eigen::VectorXd&, double)>>
Eigen::VectorXd fixedPointSolvePLaplacian(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    FUNCTOR_F f, double p, double rtol = 1.0E-5, double atol = 1.0E-8,
    unsigned int itmax = 1000,
    RECORDER&& rec = [](const Eigen::VectorXd&, double) -> void {}) {
  // Wrap right hand side source function into a Mesh Function
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Dimension of finite element space`
  const std::size_t N_dofs(dofh.NumDofs());

  // Coefficient vector for approximate solution
  Eigen::VectorXd mu_vec(N_dofs);
  // Initial guess needs to have non-zero gradient for the first
  // linear system in the fixed-point iteration to be well-posed.
  for (int i = 0; i < N_dofs; ++i) {
    const lf::mesh::Entity& node{dofh.Entity(i)};

    const Eigen::Vector2d node_pos{
        lf::geometry::Corners(*node.Geometry()).col(0)};
    mu_vec(i) = sin(3. * M_PI * node_pos[0]) * sin(3. * M_PI * node_pos[1]);
  }

  // Main fixed-point iteration loop
  Eigen::VectorXd mu_old(N_dofs);
  double diff_norm;  // Norm of correction
  unsigned int steps = 0;
  do {
    mu_old = mu_vec;
    // Solve linear fixed-point system and return solution in mu\_vec
    fixedPointNextIt(fes_p, f, p, mu_vec);
    mu_vec = .75 * mu_vec + .25 * mu_old;
    diff_norm = (mu_vec - mu_old).norm();
    rec(mu_vec, diff_norm);
    steps++;
    //  Correction-based termination
  } while ((diff_norm >= rtol * mu_vec.norm()) && (diff_norm >= atol) &&
           (steps < itmax));
  return mu_vec;
}

/** @brief Newton-based solver for Dirchlet BVP for p-Laplacian
 *
 * @param fes_p lowest-order Lagrangian FE space (Also contains information
 * about the mesh and the local-to-global index mapping)
 * @param f right-hanmd side source function
 * @param p power of Laplacian
 * @param rtol relative error tolerance
 * @param atol absolute error tolerance
 * @param itmax max Newton iterations
 * @param rec recorder to collect data
 * @return tent function basis expansion of computed solution
 */

template <typename FUNCTOR_F, typename RECORDER = std::function<
                                  void(const Eigen::VectorXd&, double)>>
Eigen::VectorXd newtonSolvePLaplacian(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    FUNCTOR_F f, double p, double rtol = 1.0E-5, double atol = 1.0E-10,
    unsigned int itmax = 100,
    RECORDER&& rec = [](const Eigen::VectorXd&, double) -> void {}) {
  // Wrap right hand side source function into a Mesh Function
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Dimension of finite element space`
  const std::size_t N_dofs(dofh.NumDofs());

  // Coefficient vector for approximate solution
  Eigen::VectorXd mu_vec(N_dofs);
  // Initial guess needs to have non-zero gradient for the first
  // linear system in the Newton iteration to be well-posed.
  for (int i = 0; i < N_dofs; ++i) {
    const lf::mesh::Entity& node{dofh.Entity(i)};

    const Eigen::Vector2d node_pos{
        lf::geometry::Corners(*node.Geometry()).col(0)};
    mu_vec(i) = sin(3. * M_PI * node_pos[0]) * sin(3. * M_PI * node_pos[1]);
  }
  // Main fixed-point iteration loop
  Eigen::VectorXd mu_old(N_dofs);
  double diff_norm;  // Norm of correction
  unsigned int steps = 0;
  do {
    mu_old = mu_vec;
    // Solve linear fixed-point system and return solution in mu\_vec
    mu_vec += newtonUpdate(fes_p, f, p, mu_vec);
    diff_norm = (mu_vec - mu_old).norm();
    rec(mu_vec, diff_norm);
    steps++;
    //  Correction-based termination
  } while ((diff_norm >= rtol * mu_vec.norm()) && (diff_norm >= atol) &&
           (steps < itmax));
  return mu_vec;
}

}  // namespace plaplacian
