#ifndef CONVECTION_EMP_H
#define CONVECTION_EMP_H

/**
 * @file convection_emp.h
 * @brief EMP for a convection term based on linear FE and the trapezoidal rule
 * @author Philippe Peter
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>

namespace ConvectionDiffusion {

/**
 * @headerfile convection_emp.h
 * @brief Computes the local matrices for the convection term based on linear
 * finite elements and the trapezoidal rule (standard Galerkin approach).
 *
 * @tparam FUNCTOR function that defines the vector valued velocity
 * coefficient v.
 */
template <typename FUNCTOR>
class ConvectionElementMatrixProvider {
 public:
  /**
   * @brief
   * @param v functor for the velocity field
   */
  explicit ConvectionElementMatrixProvider(FUNCTOR v) : v_(v) {}

  /**
   * @brief main routine for the computation of element matrices
   * @param enttity reference to the TRIANGULAR cell for which the element
   * matrix should be computed.
   * @return a 3x3 matrix containing the element matrix.
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);

  /** @brief Default implementation: all cells are active */
  bool isActive(const lf::mesh::Entity& /*entity*/) const { return true; }

 private:
  FUNCTOR v_;  // functor for the velocity field.
};

template <typename FUNCTOR>
Eigen::Matrix3d ConvectionElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);

  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  // evaluate velocity field at the vertices:
  Eigen::MatrixXd velocities(2, 3);
  velocities << v_(corners.col(0)), v_(corners.col(1)), v_(corners.col(2));

  // compute local matrix using local trapezoidal rule:
  loc_mat = lf::geometry::Volume(*geo_ptr) / 3.0 * velocities.transpose() *
            grad_basis;

  return loc_mat;
}

}  // namespace ConvectionDiffusion

#endif  // CONVECTIOIN_EMP_H