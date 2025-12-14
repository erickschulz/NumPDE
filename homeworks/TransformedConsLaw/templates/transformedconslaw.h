/**
 * @ file transformedconslaw.h
 * @ brief NPDE homework about conservation law with non-linear density
 * @ author Ralf Hiptmair, Oliver Rietmann
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <utility>

namespace TRFCL {

/**
 * @brief Class describing a Cauchy problem with non-linear density
 */
class NonStdCauchyProblemCL {
 public:
  NonStdCauchyProblemCL() = default;
  virtual ~NonStdCauchyProblemCL() = default;
  // Function rho
  double rho(double z) const;
  // Derivative of function rho
  double drho(double z) const;
  // Function g
  double g(double z) const;
  // Derivative of function g
  double dg(double z) const;
  // Finite interval containing "interesting" parts of solution
  std::pair<double, double> domain() const;
  // Final time
  double T() const;
  // Initial data
  double z0(double x) const;
};

/**
 * @brief Newton's method for the inversion of rho
 *
 * Approximately solves rho(z) = u employing Newton's method with initial guess.
 *
 * @tparam RHOFUNCTOR provides rho(z), std::function<double(double)>
 * @tparam DRHOFUNCTOR provides rho'(z), std::function<double(double)>
 * @param u value for which the the inverse of rho should be computed
 * @param z0 initial guess for Newton's method
 * @param rho functor for z -> rho(z)
 * @param rhod functor providing the derivative of rho
 * @param atol absolute tolerance
 * @param rtol relative tolerance

 WRONG IMPLEMENTATION OF TEARMINATION CRITERION.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename RHOFUNCTOR, typename DRHOFUNCTOR>
double rhoInverse(double u, double z0, RHOFUNCTOR&& rho, DRHOFUNCTOR&& drho,
                  double atol = 1.0E-10, double rtol = 1.0E-5) {
  //====================
  // Your code goes here
  //====================
  return z0;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <class CAUCHYPROBLEM>
Eigen::VectorXd semiDiscreteRhs(const Eigen::VectorXd& mu,
                                const Eigen::VectorXd& zeta,
                                CAUCHYPROBLEM prb) {
  int N = mu.size();
  Eigen::VectorXd rhs(N);

  //====================
  // Your code goes here
  //====================

  return rhs;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <class CAUCHYPROBLEM, typename RECORDER = std::function<
                                   void(double, const Eigen::VectorXd&)>>
Eigen::VectorXd solveCauchyPrb(
    unsigned int M, unsigned int N, CAUCHYPROBLEM prb,
    RECORDER&& rec = [](double /*time*/, const Eigen::VectorXd&
                        /*zstate*/) -> void {}) {
  // Get inital data for zeta
  std::pair<double, double> limits = prb.domain();
  Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, limits.first, limits.second);
  auto z0 = [&prb](double y) { return prb.z0(y); };
  Eigen::VectorXd zeta = x.unaryExpr(z0);

  // Compute time step
  double T = prb.T();
  double dt = T / M;

  // Wrap the involved functions for simpler use below
  auto rho = [&prb](double v) { return prb.rho(v); };
  auto drho = [&prb](double v) { return prb.drho(v); };
  auto r = [rho, drho](double u, double z) {
    return rhoInverse(u, z, rho, drho);
  };

  // Get inital data for mu and get the time grid
  Eigen::VectorXd mu = zeta.unaryExpr(rho);
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(M + 1, 0.0, T);

  // record solution at initial time
  rec(t[0], zeta);

  for (int k = 0; k < M; ++k) {
    //====================
    // Your code goes here
    //====================

    // record solution after current time step
    rec(t[k + 1], zeta);
  }

  return zeta;
}
/* SAM_LISTING_END_3 */

}  // namespace TRFCL
