/**
 * @file upwindfinitevolume.cc
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include "upwindfinitevolume.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>

namespace UpwindFiniteVolume {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3>& triangle) {
  Eigen::Matrix3d X;
  // Solve for the coefficients of the barycentric coordinate functions
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double computeUpwindFlux(double mui, double muk, double vhat, double dik,
                         double epsilon) {
  double flux = 0;
  //====================
  // Your code goes here
  //====================
  return flux;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::Vector2d computeCircumcenters(const Eigen::Vector2d& a1,
                                     const Eigen::Vector2d& a2,
                                     const Eigen::Vector2d& a3) {
  //====================
  // Your code goes here
  //====================
  return Eigen::Vector2d::Zero();
}
/* SAM_LISTING_BEGIN_3 */

}  // namespace UpwindFiniteVolume
