/**
 * @file systemode.h
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace SystemODE {

// Single step of RK4 for the ODE y' = f(y)
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
Eigen::VectorXd rk4step(Function&& f, double h, Eigen::VectorXd& y0) {
  Eigen::VectorXd eval(y0.size());
  //====================
  // Your code goes here
  //====================
  return eval;
}
/* SAM_LISTING_END_1 */

}  // namespace SystemODE
