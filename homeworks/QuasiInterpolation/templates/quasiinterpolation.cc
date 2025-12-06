/**
 * @file quasiinterpolation.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 15.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "quasiinterpolation.h"

#include <lf/base/base.h>  // std::span
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <memory>
#include <utility>

namespace QuasiInterpolation {

// Auxiliary function: computing the length of an edge
double edgeLength(const lf::mesh::Entity& edge) {
  Eigen::Matrix2d corners = lf::geometry::Corners(*(edge.Geometry()));
  return (corners.col(1) - corners.col(0)).norm();
}

// Auxiliary function: computing the length of the longest edge
double maxLength(const std::span<const lf::mesh::Entity* const>& edges) {
  double length = 0.0;
  for (const lf::mesh::Entity* edge : edges) {
    length = std::max(length, edgeLength(*edge));
  }
  return length;
}

/* SAM_LISTING_BEGIN_1 */
lf::mesh::utils::CodimMeshDataSet<
    std::pair<const lf::mesh::Entity*, unsigned int>>
findKp(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Variable for returning result
  lf::mesh::utils::CodimMeshDataSet<
      std::pair<const lf::mesh::Entity*, unsigned int>>
      KpMeshDataSet(mesh_p, 2);
  // Auxiliary array storing size of largest triangle adjacent to a node
  lf::mesh::utils::CodimMeshDataSet<double> sizeMeshDataSet(mesh_p, 2);
  // loop over all cells
  for (const lf::mesh::Entity* triangle : mesh_p->Entities(0)) {
    LF_ASSERT_MSG(triangle->RefEl() == lf::base::RefEl::kTria(),
                  "Only implemented for triangles");
    //====================
    // Your code goes here
    //====================
  }  // end of loop over triangles
  return KpMeshDataSet;
}
/* SAM_LISTING_END_1 */

}  // namespace QuasiInterpolation
