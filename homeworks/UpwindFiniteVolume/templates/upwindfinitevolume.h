#ifndef UPWINDFINITEVOLUME_H
#define UPWINDFINITEVOLUME_H
/**
 * @file upwindfinitevolume.h
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <cmath>

namespace UpwindFiniteVolume {

/**
 * @brief Get the coefficients of the barycentric coordinate functions for
 * a TRIA element.
 *
 * @param triangle Corners of the element.
 * @return Matrix providing the coefficients.
 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3>& triangle);

/**
 * @brief Compute the upwind flux $J_{ik}(\mu_i, \mu_k)$
 *
 * @param mui Value of $u_N$ at $p_i$
 * @param muk Value of $u_N$ at $p_k$
 * @param vhat Length of projection of $v$ onto $p_k - p_i$
 * @param dik Distance between $p_i$ and $p_k$
 * @param epsilon Strength of the diffusion
 */
double computeUpwindFlux(double mui, double muk, double vhat, double dik,
                         double epsilon);

/**
 * @brief Compute the circumcenter of a triangle.
 *
 * @param a1, a2, a3 Corners of the triangle.
 * @return Vector2d describing the circumcenter.
 */
Eigen::Vector2d computeCircumcenters(const Eigen::Vector2d& a1,
                                     const Eigen::Vector2d& a2,
                                     const Eigen::Vector2d& a3);

template <typename FUNCTOR>
class ElementMatrixProvider {
 public:
  explicit ElementMatrixProvider(FUNCTOR v, double eps) : v_(v), eps_(eps) {}

  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity& /*entity*/) const { return true; }

 private:
  FUNCTOR v_;
  double eps_;
};

/**
 * @brief Provider for the element matrix.
 *
 * @param entity Refenence to a triangular cell.
 * @return Matrix3d The element matrix.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
Eigen::Matrix3d ElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& entity) {
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

  //====================
  // Your code goes here
  //====================
  return A;
}
/* SAM_LISTING_END_1 */

template <typename FUNCTOR>
class ElementVectorProvider {
 public:
  explicit ElementVectorProvider(FUNCTOR f) : f_(f) {}

  Eigen::Vector3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity& /*entity*/) const { return true; }

 private:
  FUNCTOR f_;
};

/**
 * @brief Provider for the element vector.
 *
 * @param entity Refenence to a triangular cell.
 * @return Vector3d The element vector.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::Vector3d ElementVectorProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& entity) {
  //====================
  // Your code goes here
  //====================
  return Eigen::Vector3d::Zero();
}
/* SAM_LISTING_END_2 */
}  // namespace UpwindFiniteVolume

#endif  // define UPWINDFINITEVOLUME_H
