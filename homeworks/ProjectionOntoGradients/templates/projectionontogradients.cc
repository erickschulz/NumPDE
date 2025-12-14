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

  //====================
  // Your code goes here
  //====================
  return loc_mat;
}
/* SAM_LISTING_END_2 */

}  // namespace ProjectionOntoGradients
