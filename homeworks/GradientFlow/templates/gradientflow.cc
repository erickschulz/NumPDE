/**
 * @file gradientflow.h
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "gradientflow.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace GradientFlow {

/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd ButcherMatrix() {
  Eigen::MatrixXd A(6, 5);
  // clang-format off
  A <<       0.25,          0.,       0.,       0.,   0.,
    0.5,        0.25,       0.,       0.,   0.,
    17./50.,     -1./25.,     0.25,       0.,   0.,
    371./1360., -137./2720., 15./544.,     0.25,   0.,
    25./24.,    -49./48., 125./16., -85./12., 0.25,
    25./24.,    -49./48., 125./16., -85./12., 0.25;
  // clang-format on
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::VectorXd> SolveGradientFlow(const Eigen::VectorXd& d,
                                               double lambda,
                                               const Eigen::VectorXd& y0,
                                               double T, unsigned int M) {
  // initialize solution vector
  std::vector<Eigen::VectorXd> sol(M + 1, Eigen::VectorXd::Zero(y0.size()));

  //====================
  // Your code goes here
  //====================
  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace GradientFlow
