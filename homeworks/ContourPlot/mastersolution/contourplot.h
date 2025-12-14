/**
 * @file contourplot.h
 * @brief NPDE homework ContourPlot code
 * @author Oliver Rietmann
 * @date 25.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "ode45.h"

namespace ContourPlot {

Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg();

/* SAM_LISTING_BEGIN_0 */
template <typename GradientFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePoints(
    GradientFunctor&& gradF, Eigen::Vector2d y0, double T) {
  Eigen::Matrix<double, 2, Eigen::Dynamic> states;

  // Right-hand-side vector field of isoline ODE
  auto rhs = [gradF](Eigen::Vector2d x) -> Eigen::Vector2d {
    Eigen::Vector2d gradFx = gradF(x);
    return Eigen::Vector2d(-gradFx(1), gradFx(0)) / gradFx.norm();
  };

  // Adaptive explicit Runge-Kutta method
  Ode45<Eigen::Vector2d> integrator(rhs);
  // Set tolerances for timestep control (optional)
  integrator.options.rtol = 1e-6;
  integrator.options.atol = 1e-8;
  // Perform explicit timestepping
  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  // Convert output into requested format: points on isoline arranged into the
  // columns of a matrix.
  int M = sol.size() - 1;
  states = Eigen::MatrixXd::Zero(2, M + 1);
  for (int m = 0; m <= M; ++m) {
    states.col(m) = sol[m].first;
  }

  return states;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ(
    FFunctor&& F, Eigen::Vector2d y0, double T) {
  // Get gradient of F by symmetric difference quotients with
  // a span given by the root of the machine precision in order
  // to curb cancellation as much as possible.
  double h = std::sqrt(std::numeric_limits<double>::epsilon());
  auto gradF = [h, F](Eigen::Vector2d x) -> Eigen::Vector2d {
    Eigen::Vector2d dx(h, 0.0);
    Eigen::Vector2d dy(0.0, h);
    return Eigen::Vector2d(F(x + dx) - F(x - dx), F(x + dy) - F(x - dy)) /
           (2.0 * h);
  };

  return computeIsolinePoints(gradF, y0, T);
}
/* SAM_LISTING_END_1 */

}  // namespace ContourPlot
