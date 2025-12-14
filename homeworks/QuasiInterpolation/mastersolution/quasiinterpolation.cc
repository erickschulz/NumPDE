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
    // Fetch coordinates of vertices
    const Eigen::MatrixXd corners{lf::geometry::Corners(*triangle->Geometry())};
    // Determine size of triangle
    const double newSize = std::max({(corners.col(1) - corners.col(0)).norm(),
                                     (corners.col(2) - corners.col(1)).norm(),
                                     (corners.col(0) - corners.col(2)).norm()});
    // Obtain array of pointers to vertex objects of current triangle
    std::span<const lf::mesh::Entity* const> vertices{triangle->SubEntities(2)};
    // Loop over vertices and update size of largest adjacent triangle.
    for (unsigned int i = 0; i < 3; ++i) {
      // Note that 'size' is a reference!
      double& size = sizeMeshDataSet(*vertices[i]);
      // Current triangle is larger than those recorded earlier
      if (newSize > size) {
        // Update entry of auxiliary array
        size = newSize;
        // Store pointer to current triangle and local vertex index
        KpMeshDataSet(*vertices[i]) = std::make_pair(triangle, i);
      }
    }
  }  // end of loop over triangles
  return KpMeshDataSet;
}
/* SAM_LISTING_END_1 */

}  // namespace QuasiInterpolation
