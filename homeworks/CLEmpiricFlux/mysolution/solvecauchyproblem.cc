/**
 * @file solvecauchyproblem.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 19.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "solvecauchyproblem.h"

#include <Eigen/Core>
#include <cmath>

#include "uniformcubicspline.h"

namespace CLEmpiricFlux {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d findSupport(const UniformCubicSpline& f,
                            Eigen::Vector2d initsupp, double t) {
  Eigen::Vector2d result;
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveCauchyProblem(const UniformCubicSpline& f,
                                   const Eigen::VectorXd& mu0, double h,
                                   double T) {
  Eigen::VectorXd muT(mu0.size());
  //====================
  // Your code goes here
  //====================
  return muT;
}
/* SAM_LISTING_END_4 */

}  // namespace CLEmpiricFlux
