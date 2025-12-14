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
  const double delta = 1e-8;
  const double exponent = vhat / epsilon * dik;
  if (std::fabs(exponent) > delta) {
    // No risk of cancellation
    flux = vhat * (muk - mui) / (1 - std::exp(-exponent)) + vhat * mui;
  } else {
    // Use taylor approximation to prevent cancelation
    flux = epsilon * (muk - mui) / (dik * (1 - 0.5 * exponent)) + vhat * mui;
  }
  return flux;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::Vector2d computeCircumcenters(const Eigen::Vector2d& a1,
                                     const Eigen::Vector2d& a2,
                                     const Eigen::Vector2d& a3) {

  Eigen::Vector2d mp1 = 0.5 * (a1 + a2);
  Eigen::Vector2d mp2 = 0.5 * (a2 + a3);

  Eigen::Vector2d dir1 = a2 - a1;
  Eigen::Vector2d dir2 = a3 - a2;

  Eigen::Matrix2d dir;
  dir << dir1(1), dir2(1), -dir1(0), -dir2(0);

  Eigen::Vector2d p = dir.colPivHouseholderQr().solve(mp2 - mp1);

  Eigen::Vector2d center = mp1 + p(0) * dir.col(0);

  // Barycentric coordinates from
  // Christer Ericson’s ’Real-Time Collision Detection’

  Eigen::Vector2d v0 = a2 - a1, v1 = a3 - a1, v2 = center - a1;
  double d00 = v0.dot(v0);
  double d01 = v0.dot(v1);
  double d11 = v1.dot(v1);
  double d20 = v2.dot(v0);
  double d21 = v2.dot(v1);

  double denom = d00 * d11 - d01 * d01;
  double v = (d11 * d20 - d01 * d21) / denom;
  double w = (d00 * d21 - d01 * d20) / denom;
  double u = 1.0 - v - w;

  if (v <= 0 || w <= 0 || u <= 0) {
    throw std::runtime_error("Obtused triangle!");
    return Eigen::Vector2d::Zero();
  }
  return center;

}
/* SAM_LISTING_BEGIN_3 */

}  // namespace UpwindFiniteVolume
