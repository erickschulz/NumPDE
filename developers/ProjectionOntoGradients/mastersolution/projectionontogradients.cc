/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author Erick Schulz, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include "projectionontogradients.h"

namespace ProjectionOntoGradients {

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d ElementMatrixProvider::Eval(const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

#if SOLUTION
  // get area of the entity
  const double area = lf::geometry::Volume(*geo_ptr);

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  loc_mat = area * (grad_basis.transpose() * grad_basis);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return loc_mat;
}
/* SAM_LISTING_END_2 */

}  // namespace ProjectionOntoGradients
