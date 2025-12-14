#ifndef SEMIMPRK_H_
#define SEMIMPRK_H_

/**
 * @file semimprk.h
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>

namespace SemImpRK {

// Solve the autonomous IVP y' = f(y), y(0) = y0 using the Rosenbrock method
/* SAM_LISTING_BEGIN_0 */
template <class Func, class Jac>
std::vector<Eigen::VectorXd> SolveRosenbrock(Func&& f, Jac&& df,
                                             const Eigen::VectorXd& y0,
                                             unsigned int M, double T) {
  // Will contain all states computed by the ROW-SSM
  std::vector<Eigen::VectorXd> res(M + 1);
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_0 */

double CvgRosenbrock();

}  // namespace SemImpRK

#endif  // #ifndef SEMIMPRK_H_
