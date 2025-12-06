/**
 * @file solvecauchyproblem.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 19.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "solvecauchyproblem.h"

#include <Eigen/Core>
#include <cmath>

#include "uniformcubicspline.h"

namespace CLEmpiricFlux {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d findSupport(const UniformCubicSpline& f,
                            Eigen::Vector2d initsupp, double t) {
  Eigen::Vector2d result;
  Eigen::Vector2d speed = {f.derivative(-1.0), f.derivative(1.0)};
  result = initsupp + t * speed;
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd semiDiscreteRhs(const Eigen::VectorXd& mu0, double h,
                                FUNCTOR&& numFlux) {
  int m = mu0.size();
  Eigen::VectorXd mu1(m);
  mu1(0) = -1.0 / h * (numFlux(mu0(0), mu0(1)) - numFlux(mu0(0), mu0(0)));
  for (int j = 1; j < m - 1; ++j) {
    mu1(j) =
        -1.0 / h * (numFlux(mu0(j), mu0(j + 1)) - numFlux(mu0(j - 1), mu0(j)));
  }
  mu1(m - 1) =
      -1.0 / h *
      (numFlux(mu0(m - 1), mu0(m - 1)) - numFlux(mu0(m - 2), mu0(m - 1)));
  return mu1;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::VectorXd RalstonODESolver(FUNCTOR&& rhs, Eigen::VectorXd mu0, double tau,
                                 int n) {
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd k1 = rhs(mu0);
    Eigen::VectorXd k2 = rhs(mu0 + tau * 2.0 / 3.0 * k1);
    mu0 = mu0 + 0.25 * tau * (k1 + 3.0 * k2);
  }
  return mu0;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveCauchyProblem(const UniformCubicSpline& f,
                                   const Eigen::VectorXd& mu0, double h,
                                   double T) {
  Eigen::VectorXd muT(mu0.size());
  double tau = std::min(h / std::abs(f.derivative(-1.0)),
                        h / std::abs(f.derivative(1.0)));
  double n = (int)std::floor(T / tau);
  GodunovFlux godunovFlux(f);
  auto rhs = [h, &godunovFlux](const Eigen::VectorXd& mu) {
    return semiDiscreteRhs(mu, h, godunovFlux);
  };
  muT = RalstonODESolver(rhs, mu0, tau, n);
  return muT;
}
/* SAM_LISTING_END_4 */

}  // namespace CLEmpiricFlux
