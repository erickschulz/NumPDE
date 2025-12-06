#ifndef MIRK_H_
#define MIRK_H_

/**
 * @file mirk.h
 * @brief NPDE homework MIRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>

namespace MIRK {

/** Perform 2 steps of the Newton method applied to F and its Jacobian DF */
/* SAM_LISTING_BEGIN_0 */
template <class Func, class Jac>
Eigen::VectorXd Newton2Steps(Func&& F, Jac&& DF, Eigen::VectorXd z) {
  //====================
  // Your code goes here
  //====================
  return z;
}
/* SAM_LISTING_END_0 */

/** Perform a single step of the MIRK scheme applied to the scalar ODE
 * y' = f(y) */
/* SAM_LISTING_BEGIN_1 */
template <class Func, class Jac>
double MIRKStep(Func&& f, Jac&& df, double y0, double h) {
  // Coefficients of MIRK
  const double v1 = 1.0;
  const double v2 = 344.0 / 2025.0;
  const double d21 = -164.0 / 2025.0;
  const double b1 = 37.0 / 82.0;
  const double b2 = 45.0 / 82.0;

  //====================
  // Your code goes here
  // TODO: Replace the following dummy return value
  //====================
  return 0.0;
}
/* SAM_LISTING_END_1 */

/** Solve an ODE y' = f(y) using MIRK scheme on equidistant steps,
 * return the approximation of y(T)*/
/* SAM_LISTING_BEGIN_2 */
template <class Func, class Jac>
double MIRKSolve(Func&& f, Jac&& df, double y0, double T, unsigned int M) {
  //====================
  // Your code goes here
  //====================
  return 0.0;
}
/* SAM_LISTING_END_2 */

}  // namespace MIRK

#endif  // #ifndef MIRK_H_
