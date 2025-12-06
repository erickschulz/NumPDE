/**
 * @file leastsquaresadvection.cc
 * @brief NPDE homework LeastSquaresAdvection code
 * @author Ralf Hiptmair
 * @date July 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#include "leastsquaresadvection.h"

#include <Eigen/src/Core/util/Constants.h>
#include <lf/base/lf_assert.h>
#include <lf/base/types.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/special_entity_sets.h>

namespace LeastSquaresAdvection {
lf::mesh::utils::AllCodimMeshDataSet<bool> flagEntitiesOnInflow(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    Eigen::Vector2d velocity) {
  LF_ASSERT_MSG(mesh_p->DimMesh() == 2, "Only implemented for 2D meshes!");
  const lf::mesh::utils::CodimMeshDataSet<bool> bd_ed_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  lf::mesh::utils::AllCodimMeshDataSet<bool> inflow_flags{mesh_p, false};
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    const lf::base::size_type num_vert = cell->RefEl().NumNodes();
    const Eigen::Matrix<double, 2, Eigen::Dynamic> corners{
        lf::geometry::Corners(*cell->Geometry())};
    // Loop over the edges of the cell
    std::span<const lf::mesh::Entity* const> edges{cell->SubEntities(1)};
    const lf::base::size_type num_edges = cell->RefEl().NumSubEntities(1);
    LF_ASSERT_MSG(edges.size() == num_edges, "Num edges mismatch!");
    LF_ASSERT_MSG(num_edges == num_vert,
                  "No. od edges/vertices must be the same");
    for (int edge_idx = 0; edge_idx < num_edges; edge_idx++) {
      // Check whether edge is on the boundasry
      if (bd_ed_flags(*edges[edge_idx])) {
        // Obtain vector in the direction of the edge
        const int p0_idx = edge_idx;
        const int p1_idx = (edge_idx + 1) % num_vert;
        const int pext_idx = (edge_idx + 2) % num_vert;
        const Eigen::Vector2d ed_dir =
            corners.col(p1_idx) - corners.col(p0_idx);
        // Normal vector
        Eigen::Vector2d ed_n(-ed_dir[1], ed_dir[0]);
        // Ensure that the normal vector points to the exterior of the cell
        if (ed_n.dot(corners.col(pext_idx) - corners.col(p0_idx)) > 0.0) {
          ed_n = -ed_n;
        }
        // Compare direction of velocity and that of the exterior normal
        if (ed_n.dot(velocity) <= 0.0) {
          inflow_flags(*edges[edge_idx]) = true;
          for (const lf::mesh::Entity* node : edges[edge_idx]->SubEntities(1)) {
            inflow_flags(*node) = true;
          }
        }
      }
    }
  }
  return inflow_flags;
}

}  // namespace LeastSquaresAdvection
