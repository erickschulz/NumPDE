/**
 * @file semilinearellipticbvp.h
 * @brief NPDE homework SemilinearEllipticBVP code
 * @author R. Hiptmair
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

namespace semilinearellipticbvp {

/** @brief Auxiliary mesh function for wrapping a function around a mesh
 * function
 *
 * @tparam MESHFUNCTION the base mesh function, which is wrapped
 * @tparam FUNCTION function type providing the wrapping
 */
/* SAM_LISTING_BEGIN_2 */
template <typename MESHFUNCTION, typename FUNCTION>
class FunctionMFWrapper {
 public:
  static_assert(lf::mesh::utils::MeshFunction<MESHFUNCTION>);
  using mf_vector_t = decltype(std::declval<MESHFUNCTION>()(
      std::declval<const lf::mesh::Entity>(),
      std::declval<const Eigen::MatrixXd>()));
  using mf_result_t = typename mf_vector_t::value_type;
  // Return type of the function F
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

/** @brief Computation of fixed-point update
 *
 * This function computes the update for the fixed-point iteration, which is an
 * iterative solver for the semilinear elliptic BVP.
 */
void fixedPointNextIt(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    Eigen::VectorXd& mu_vec, Eigen::VectorXd& rhs_vec);

/** @brief Computation of Newton iteration update
 *
 * This function computes the update for the Newton iteration, which is an
 * iterative solver for the semilinear elliptic BVP.
 */
void newtonNextIt(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    Eigen::VectorXd& mu_vec, Eigen::VectorXd& rhs_vec);

/** @brief Solves semi-linear elliptic BVP by means of a fixed-point iteration
 *
 * @param fes_p pointer to discrete trial and test finite-element space
 * @param f functor object for right hand side source function
 * @param rtol relative tolerance for correction-based termination criterion
 * @param atol absolute tolerance for correction-based termination criterion
 */

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR_F, typename RECORDER = std::function<
                                  void(const Eigen::VectorXd&, double)>>
Eigen::VectorXd solveSemilinearBVP(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    FUNCTOR_F f, double rtol = 1.0E-4, double atol = 1.0E-8,
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

  // Assemble right-hand side vector
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  // Assemble volume part of right-hand side vector depending on the source
  // function f.
  // Initialize object taking care of local computations on all cells.
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fes_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Coefficient vector for approximate solution
  Eigen::VectorXd mu_vec(N_dofs);
  // Initial guesd = 0
  mu_vec.setZero();
  // Main fixed-point iteration loop
  Eigen::VectorXd mu_old(N_dofs);
  double diff_norm;  // Norm of correction
  unsigned int steps = 0;
  do {
    mu_old = mu_vec;
    // Solve linear fixed-point system and return solution in mu\_vec
    fixedPointNextIt(fes_p, mu_vec, phi);
    diff_norm = (mu_vec - mu_old).norm();
    rec(mu_vec, diff_norm);
    steps++;
    // Correction-based termination
  } while ((diff_norm >= rtol * mu_vec.norm()) && (diff_norm >= atol) &&
           (steps < itmax));
  return mu_vec;
}
/* SAM_LISTING_END_1 */

/** @brief Solves semi-linear elliptic BVP by means of a Newton iteration
 *
 * @param fes_p pointer to discrete trial and test finite-element space
 * @param f functor object for right hand side source function
 * @param rtol relative tolerance for correction-based termination criterion
 * @param atol absolute tolerance for correction-based termination criterion
 */

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR_F, typename RECORDER = std::function<
                                  void(const Eigen::VectorXd&, double)>>
Eigen::VectorXd newtonSolveSemilinearBVP(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    FUNCTOR_F f, double rtol = 1.0E-4, double atol = 1.0E-8,
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

  // Assemble right-hand side vector
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  // Assemble volume part of right-hand side vector depending on the source
  // function f.
  // Initialize object taking care of local computations on all cells.
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fes_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Coefficient vector for approximate solution
  Eigen::VectorXd mu_vec(N_dofs);
  // Initial guesd = 0
  mu_vec.setZero();
  // Main Newton iteration loop
  Eigen::VectorXd mu_old(N_dofs);
  double diff_norm;  // Norm of correction
  unsigned int steps = 0;
  do {
    mu_old = mu_vec;
    // Solve linear Newton system and return solution in mu\_vec
    newtonNextIt(fes_p, mu_vec, phi);
    diff_norm = (mu_vec - mu_old).norm();
    rec(mu_vec, diff_norm);
    steps++;
    // Correction-based termination
  } while ((diff_norm >= rtol * mu_vec.norm()) && (diff_norm >= atol) &&
           (steps < itmax));
  return mu_vec;
}
/* SAM_LISTING_END_2 */

void testSolverSemilinearBVP(unsigned int reflevels = 3);

}  // namespace semilinearellipticbvp
