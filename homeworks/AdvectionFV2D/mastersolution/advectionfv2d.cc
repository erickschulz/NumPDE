/**
 * @file advectionfv2d.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include "advectionfv2d.h"

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace AdvectionFV2D {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3>& triangle) {
  Eigen::Matrix3d X;

  // solve for the coefficients of the barycentric coordinate functions
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>>>
computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Appears only in mastersolution

  // Initialize datastruture for the result
  lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>>
      result(mesh_p, 0);

  // Compute normal vectors
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    if (corners.cols() == 3) {
      Eigen::Matrix<double, 2, Eigen::Dynamic> normal_vectors(2, 3);

      // Compute normal vectors
      Eigen::Matrix<double, 2, 3> bary_coord = gradbarycoordinates(corners);

      // Reorder computed normal vectors
      // normal_vectors[0] -> normal_vector of edge[0]
      normal_vectors.col(0) = -(bary_coord.col(2)).normalized();
      normal_vectors.col(1) = -(bary_coord.col(0)).normalized();
      normal_vectors.col(2) = -(bary_coord.col(1)).normalized();

      // Save the matrix in our datastructure
      result(*cell) = normal_vectors;
    } else {  // corners.cols() == 4
      Eigen::Matrix<double, 2, Eigen::Dynamic> normal_vectors(2, 4);

      // Split the quadrilateral into two triangles (1,2,3) and (1,3,4)
      Eigen::Matrix<double, 2, 3> tria_1;
      Eigen::Matrix<double, 2, 3> tria_2;
      tria_1 << corners.col(0), corners.col(1), corners.col(2);
      tria_2 << corners.col(0), corners.col(2), corners.col(3);

      // Compute normal vectors
      Eigen::Matrix<double, 2, 3> bary_coord_1 = gradbarycoordinates(tria_1);
      Eigen::Matrix<double, 2, 3> bary_coord_2 = gradbarycoordinates(tria_2);

      // Reorder computed normal vectors
      // normal_vectors[0] -> normal_vector of edge[0]
      normal_vectors.col(0) = -(bary_coord_1.col(2)).normalized();
      normal_vectors.col(1) = -(bary_coord_1.col(0)).normalized();
      normal_vectors.col(2) = -(bary_coord_2.col(0)).normalized();
      normal_vectors.col(3) = -(bary_coord_2.col(1)).normalized();

      // Save the matrix in our datastructure
      result(*cell) = normal_vectors;
    }
  }
  return std::make_shared<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>(result);

}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
getAdjacentCellPointers(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Initialize auxilary object
  lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 2>>
      aux_obj(mesh_p, 1, {nullptr, nullptr});

  // Iterate over every cell
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    auto cell_edges = cell->SubEntities(1);

    // Iterate over every edge of the cell
    for (const lf::mesh::Entity* edge : cell_edges) {
      // If aux_obj at index 0 was not set, save cell there
      // otherwise save at the second position
      // (The first position has to be set; the second might be set)
      if (aux_obj(*edge)[0] == nullptr) {
        aux_obj(*edge)[0] = cell;
      } else if (aux_obj(*edge)[1] == nullptr) {
        aux_obj(*edge)[1] = cell;
      } else {
        throw std::runtime_error("Error in aux_obj");
      }
    }
  }

  // Initialize datastructure for result
  lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>
      result(mesh_p, 0, {nullptr, nullptr, nullptr, nullptr});

  // Collect the objects of the auxilary object
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    auto cell_edges = cell->SubEntities(1);
    int counter = 0;
    for (const lf::mesh::Entity* edge : cell_edges) {
      if (aux_obj(*edge)[0] != cell) {
        result(*cell)[counter] = aux_obj(*edge)[0];
      } else {
        result(*cell)[counter] = aux_obj(*edge)[1];
      }
      counter++;
    }
  }

  return std::make_shared<lf::mesh::utils::CodimMeshDataSet<
      std::array<const lf::mesh::Entity*, 4>>>(result);
}
/* SAM_LISTING_END_3 */

// Function returning the barycenter of TRIA or QUAD
/* SAM_LISTING_BEGIN_4 */
Eigen::Vector2d barycenter(const Eigen::MatrixXd corners) {
  Eigen::Vector2d midpoint;
  if (corners.cols() == 3) {
    midpoint = (corners.col(0) + corners.col(1) + corners.col(2)) / 3.0;
  } else if (corners.cols() == 4) {
    midpoint =
        (corners.col(0) + corners.col(1) + corners.col(2) + corners.col(3)) /
        4.0;
  } else {
    throw std::runtime_error("Wrong geometrie in barycenter()");
  }
  return midpoint;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
double computeHmin(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Vector to store the distances between cells
  std::vector<double> min_h;

  // Get Adjectent Cells
  std::shared_ptr<
      lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  // Iterate over all cells
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Compute the barycenter of the cell
    Eigen::Vector2d cur_midpoint = barycenter(corners);

    // Iterate over all adjecent cells
    for (const lf::mesh::Entity* neighbour_cell : (*adjacentCells)(*cell)) {
      // Check that the neighbor exists
      if (neighbour_cell != nullptr) {
        const lf::geometry::Geometry* geo_p_neighbour =
            neighbour_cell->Geometry();
        const Eigen::MatrixXd neighbour_corners =
            lf::geometry::Corners(*geo_p_neighbour);

        // Compute barycenter of neighbour cell
        Eigen::Vector2d neighbour_midpoint = barycenter(neighbour_corners);

        // Compute distances between cells
        double distance = (cur_midpoint - neighbour_midpoint).norm();

        // Store value in vector
        min_h.push_back(distance);
      }
    }
  }

  // Find the minimum value in the vector and return it
  return *std::min_element(std::begin(min_h), std::end(min_h));
}
/* SAM_LISTING_END_5 */

}  // namespace AdvectionFV2D
