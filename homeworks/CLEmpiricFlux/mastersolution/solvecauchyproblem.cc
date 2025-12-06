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
