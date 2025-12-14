#ifndef CROSSPROD_H_
#define CROSSPROD_H_

/**
 * @file crossprod.h
 * @brief NPDE homework CrossProd code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "implicitrkintegrator.h"

namespace CrossProd {

/* SAM_LISTING_BEGIN_1 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_imp_mid(Function&& f, Jacobian&& Jf,
                                           double T, const Eigen::VectorXd& y0,
                                           unsigned int M) {
  std::vector<Eigen::VectorXd> res(M + 1);
  // Construct the implicit mid-point method with the class
  // implicit_RKIntegrator and execute the .solve() method.
  // Return the vector containing all steps including initial and final value.
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_lin_mid(Function&& f, Jacobian&& Jf,
                                           double T, const Eigen::VectorXd& y0,
                                           unsigned int M) {
  std::vector<Eigen::VectorXd> res;
  // Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing
  // all steps including initial and final value.
  //====================
  // Your code goes here
  //====================

  return res;
}
/* SAM_LISTING_END_2 */

void tab_crossprod();

}  // namespace CrossProd

#endif  // #define CROSSPROD
