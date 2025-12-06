/** @file
 * @brief Extra functions for homework problems StationaryCurrents
 * @author Ralf Hiptmair
 * @date July 2020
 * @copyright MIT License
 */

#include "stationarycurrents_supplement.h"

#include "stationarycurrents.h"

namespace dmxbc {
// Debugging function
void printNodeTags(const lf::mesh::Mesh& mesh,
                   lf::mesh::utils::CodimMeshDataSet<int>& nodeids) {
  // Container for counters
  std::map<int, unsigned int> counters;
  // Loop over nodes of the mesh and counter occurrence of ids
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    if (nodeids.DefinedOn(*node)) {
      counters[nodeids(*node)]++;
    } else {
      std::cout << " Node " << *node << " has no id!" << std::endl;
    }
  }
  for (auto& cnt : counters) {
    std::cout << "id = " << cnt.first << ": " << cnt.second << " nodes"
              << std::endl;
  }
}  // end printNodeTags

std::tuple<Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 2, 3>, double>
getTriangleGradLambdaNormals(Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9 in the lecture document
  Eigen::Matrix<double, 3, 3> X;
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix grad\_bary\_coords
  const auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};
  // The determinant of the auxiliary matrix also supplies the area of the
  // triangle
  const double twicearea = std::abs(X.determinant());
  return {grad_bary_coords,
          -twicearea * grad_bary_coords *
              Eigen::PermutationMatrix<3>(Eigen::Vector3i(2, 0, 1)),
          twicearea / 2};
}  // end getTriangleGradLambdaNormals

Eigen::Matrix<double, 2, 3> exteriorTriangleNormals(
    Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9
  Eigen::Matrix<double, 3, 3> X;
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix grad\_bary\_coords
  const auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};
  // The determinant of the auxiliary matrix also supplies the area of the
  // triangle
  const double twicearea = std::abs(X.determinant());
  return -twicearea * grad_bary_coords *
         Eigen::PermutationMatrix<3>(Eigen::Vector3i(2, 0, 1));
}  // end exteriorTriangleNormals

Eigen::MatrixXd exteriorCellNormals(const Eigen::MatrixXd& corners) {
  LF_ASSERT_MSG(corners.rows() == 2,
                "Columns must contains coordinates of 2D points");
  // Number of vertices
  unsigned int m = corners.cols();
  // Center of gravity
  const Eigen::Vector2d c{corners.rowwise().sum() / m};
  // Matrix in whose columns we return the normals
  Eigen::MatrixXd normals(corners.rows(), corners.cols());
  for (int j = 0; j < m; ++j) {
    // Edge vector
    const Eigen::Vector2d edge_vec{corners.col((j + 1) % m) - corners.col(j)};
    // Vector from barycenter to endpoint of the edge
    const Eigen::Vector2d a{corners.col(j) - c};
    // Normal vector
    normals.col(j) = Eigen::Vector2d(-edge_vec[1], edge_vec[0]);
    // Adjust direction
    if (a.dot(normals.col(j)) < 0) {
      normals.col(j) *= -1;
    }
  }  // End: loop over straight edges
  return normals;
}  // end: exteriorCellNormals

lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> exteriorEdgeWeightedNormals(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Array indexed by edges for returning exterior normals
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> extnormals(
      mesh_p, 1, Eigen::Vector2d(0.0, 0.0));
  // Find edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    // Make sure the cell is of triangular shape
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "Only implemented for triangular cells");
    const auto vertices{lf::geometry::Corners(*(cell->Geometry()))};
    const Eigen::Matrix<double, 2, 3> normals{
        exteriorTriangleNormals(vertices)};
    // Obtain array of edge pointers (relative co-dimension = 1)
    std::span<const lf::mesh::Entity* const> sub_ent_range{
        cell->SubEntities(1)};
    // loop over the edges and check whether they belong to the boundary
    for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
      const lf::mesh::Entity& edge{*sub_ent_range[j]};
      if (bd_flags(edge)) {
        // Found edge on the boundary. Set normal vector
        extnormals(edge) = normals.col(j);
      }
    }
  }
  return extnormals;
}  // end exteriorEdgeWeightedNormals

// A debugging function
bool validateNormals(const lf::mesh::Mesh& mesh) {
  // Run through cells and compute the normals
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // NOLINTBEGIN(clang-analyzer-deadcode.DeadStores)
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    // NOLINTEND(clang-analyzer-deadcode.DeadStores)
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "implemented for triangles only");
    const lf::geometry::Geometry& geo{*(cell->Geometry())};
    const auto vertices{lf::geometry::Corners(geo)};

    auto etn{exteriorTriangleNormals(vertices)};
    auto ecn{exteriorCellNormals(vertices)};
    if ((etn - ecn).norm() > 1E-8) {
      std::cout << "Mismatch of exterior normals for " << std::endl
                << vertices << std::endl;
      std::cout << "etn = \n" << etn << std::endl;
      std::cout << "ecn = \n" << ecn << std::endl;
      return false;
    }
  }  // end loop over edges
  std::cout << "Normals ok" << std::endl;
  return true;
}  // end validate normals

}  // namespace dmxbc
