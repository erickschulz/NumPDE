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
  Eigen::VectorXd k1 = f(y0);
  Eigen::VectorXd k2 = f(y0 + h / 2 * k1);
  Eigen::VectorXd k3 = f(y0 + h / 2 * k2);
  Eigen::VectorXd k4 = f(y0 + h * k3);
  eval = y0 + h / 6 * k1 + h / 3 * k2 + h / 3 * k3 + h / 6 * k4;
  return eval;
}
/* SAM_LISTING_END_1 */

}  // namespace SystemODE
