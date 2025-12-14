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
#if SOLUTION
  // First Newton step
  z = z - DF(z).lu().solve(F(z));
  // Second Newton step
  z = z - DF(z).lu().solve(F(z));
#else
  //====================
  // Your code goes here
  //====================
#endif
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

#if SOLUTION
  // F derived from MIRK scheme (vector valued)
  auto F = [&f, &y0, &v1, &v2, &d21, &b1, &b2,
            &h](Eigen::Vector3d z) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << z(0) - (1 - v1) * y0 - v1 * z(2),
        z(1) - (1 - v2) * y0 - v2 * z(2) - h * d21 * f(z(0)),
        z(2) - y0 - h * (b1 * f(z(0)) + b2 * f(z(1)));
    return ret;
  };
  // Jacobian of F (matrix-valued)
  auto DF = [&df, &v1, &v2, &d21, &b1, &b2,
             &h](Eigen::Vector3d z) -> Eigen::Matrix3d {
    Eigen::Matrix3d M;
    M << 0, 0, v1, h * d21 * df(z(0)), 0, v2, h * b1 * df(z(0)),
        h * b2 * df(z(1)), 0;
    return Eigen::Matrix3d::Identity() - M;
  };
  // Initial Guess for Newton steps
  Eigen::Vector3d z(0.0, 0.0, y0);

  // Approximate z, s.t. F(z) = 0
  z = Newton2Steps(F, DF, z);
  return z(2);
#else
  //====================
  // Your code goes here
  // TODO: Replace the following dummy return value
  //====================
  return 0.0;
#endif
}
/* SAM_LISTING_END_1 */

/** Solve an ODE y' = f(y) using MIRK scheme on equidistant steps,
 * return the approximation of y(T)*/
/* SAM_LISTING_BEGIN_2 */
template <class Func, class Jac>
double MIRKSolve(Func&& f, Jac&& df, double y0, double T, unsigned int M) {
#if SOLUTION
  // Step size
  const double h = T / M;
  // Will contain next step
  double y = y0;
  // Perform N quidistant steps
  for (unsigned int i = 0; i < M; ++i) {
    y = MIRKStep(f, df, y, h);
  }
  // Return final value at t = T
  return y;
#else
  //====================
  // Your code goes here
  //====================
  return 0.0;
#endif
}
/* SAM_LISTING_END_2 */

}  // namespace MIRK

#endif  // #ifndef MIRK_H_
