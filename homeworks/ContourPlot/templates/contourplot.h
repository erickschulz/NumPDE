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

  //====================
  // Your code goes here
  // Use ode45 to solve the ODE with right-hand side
  // given by the function rhs(...) and fill the
  // matrix states with the solution vectors
  //====================

  return states;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ(
    FFunctor&& F, Eigen::Vector2d y0, double T) {
  //====================
  // Your code goes here
  // Replace the following dummy return value
  // by the matrix containing the isoline points:
  return Eigen::Matrix<double, 2, 42>::Zero();
  //====================
}
/* SAM_LISTING_END_1 */

}  // namespace ContourPlot
