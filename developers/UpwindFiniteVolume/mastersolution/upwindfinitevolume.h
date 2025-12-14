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

#if SOLUTION
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type!");

  const lf::geometry::Geometry* geo_p = entity.Geometry();
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

  Eigen::Vector2d circumcenter =
      computeCircumcenters(corners.col(0), corners.col(1), corners.col(2));

  Eigen::Matrix3d f = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d d = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d vv = Eigen::Matrix3d::Zero();

  f(0, 1) = (0.5 * (corners.col(0) + corners.col(1)) - circumcenter).norm();
  f(0, 2) = (0.5 * (corners.col(0) + corners.col(2)) - circumcenter).norm();
  f(1, 2) = (0.5 * (corners.col(1) + corners.col(2)) - circumcenter).norm();
  f.triangularView<Eigen::Lower>() = f.transpose();

  d(0, 1) = (corners.col(1) - corners.col(0)).norm();
  d(0, 2) = (corners.col(2) - corners.col(0)).norm();
  d(1, 2) = (corners.col(2) - corners.col(1)).norm();
  d.triangularView<Eigen::Lower>() = d.transpose();

  vv(0, 1) = v_(0.5 * (corners.col(0) + corners.col(1)))
                 .dot((corners.col(1) - corners.col(0)).normalized());
  vv(0, 2) = v_(0.5 * (corners.col(0) + corners.col(2)))
                 .dot((corners.col(2) - corners.col(0)).normalized());
  vv(1, 2) = v_(0.5 * (corners.col(1) + corners.col(2)))
                 .dot((corners.col(2) - corners.col(1)).normalized());
  vv.triangularView<Eigen::Lower>() = -vv.transpose();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i != j) {
        if (std::abs(vv(i, j)) > 1e-14) {
          double beta = 1. / (1. - std::exp(-vv(i, j) * d(i, j) / eps_));
          A(i, j) = vv(i, j) * beta * f(i, j);
          A(i, i) += vv(i, j) * (1. - beta) * f(i, j);
        } else {
          A(i, j) = eps_ / d(i, j) * f(i, j);
          A(i, i) -= eps_ / d(i, j) * f(i, j);
        }
      }
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
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
#if SOLUTION
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type!");

  const lf::geometry::Geometry* geo_p = entity.Geometry();
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

  Eigen::Vector2d circumcenter =
      computeCircumcenters(corners.col(0), corners.col(1), corners.col(2));

  double a12 = (corners.col(0) - corners.col(1)).norm() *
               (0.5 * (corners.col(0) + corners.col(1)) - circumcenter).norm() /
               4.0;
  double a13 = (corners.col(0) - corners.col(2)).norm() *
               (0.5 * (corners.col(0) + corners.col(2)) - circumcenter).norm() /
               4.0;
  double a23 = (corners.col(1) - corners.col(2)).norm() *
               (0.5 * (corners.col(1) + corners.col(2)) - circumcenter).norm() /
               4.0;

  Eigen::Vector3d result;
  result(0) = (a12 + a13) * f_(corners.col(0));
  result(1) = (a12 + a23) * f_(corners.col(1));
  result(2) = (a13 + a23) * f_(corners.col(2));

  return result;
#else
  //====================
  // Your code goes here
  //====================
  return Eigen::Vector3d::Zero();
#endif
}
/* SAM_LISTING_END_2 */
}  // namespace UpwindFiniteVolume

#endif  // define UPWINDFINITEVOLUME_H
