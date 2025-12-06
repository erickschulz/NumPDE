#ifndef EXPONENTIALINTEGRATOR_H_
#define EXPONENTIALINTEGRATOR_H_

/**
 * @file exponentialintegrator.h
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace ExponentialIntegrator {

Eigen::MatrixXd phim(const Eigen::MatrixXd& Z);

// Calculate a single step of the exponential Euler method.
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
Eigen::VectorXd exponentialEulerStep(const Eigen::VectorXd& y0, Function&& f,
                                     Jacobian&& df, double h) {
  //====================
  // Your code goes here
  // Replace the following dummy return value:
  return y0;
  //====================
}
/* SAM_LISTING_END_0 */

void testExpEulerLogODE();

}  // namespace ExponentialIntegrator

#endif  // #ifndef EXPONENTIALINTEGRATOR_H_
