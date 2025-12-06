/**
 * @file mesh.cc
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mesh.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <Eigen/Core>
#include <array>
#include <memory>

std::shared_ptr<lf::mesh::Mesh> Generate2DTestMesh() {
  using size_type = lf::mesh::Mesh::size_type;

  // Obtain mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  double scale = 1.0;

  // Hybrid mesh of square [0,3]^2 with only triangles and rectangles
  // Set coordinates of nodes
  // clang-format off
  std::array<std::array<double, 2>, 10> node_coord{
      std::array<double, 2>({0, 0 }),
      std::array<double, 2>({1, 0 }),
      std::array<double, 2>({2, 0 }),
      std::array<double, 2>({3, 0 }),
      std::array<double, 2>({0 ,2 }),
      std::array<double, 2>({1 ,2 }),
      std::array<double, 2>({3 ,2 }),
      std::array<double, 2>({0 ,3 }),
      std::array<double, 2>({1 ,3 }),
      std::array<double, 2>({3 ,3 })};
  // clang-format on

  // Specify triangles (five)
  std::array<std::array<size_type, 3>, 5> tria_nodes{
      std::array<size_type, 3>({1, 2, 5}), std::array<size_type, 3>({2, 3, 6}),
      std::array<size_type, 3>({5, 2, 6}), std::array<size_type, 3>({4, 5, 7}),
      std::array<size_type, 3>({5, 8, 7})};

  // Specify parallelograms (two)
  std::array<std::array<size_type, 4>, 2> parg_nodes{
      std::array<size_type, 4>({0, 1, 5, 4}),
      std::array<size_type, 4>({5, 6, 9, 8})};

  // Create nodes
  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(
        Eigen::Vector2d({node[0] * scale, node[1] * scale}));
  }

  // generate triangles
  for (const auto& node : tria_nodes) {
    mesh_factory_ptr->AddEntity(
        lf::base::RefEl::kTria(),
        std::span<const size_type>({node[0], node[1], node[2]}),
        std::unique_ptr<lf::geometry::Geometry>(nullptr));
  }

  // generate Parallelograms
  for (const auto& node : parg_nodes) {
    Eigen::MatrixXd quad_coord(2, 4);
    for (int n_pt = 0; n_pt < 4; ++n_pt) {
      quad_coord(0, n_pt) = node_coord[node[n_pt]][0];
      quad_coord(1, n_pt) = node_coord[node[n_pt]][1];
    }
    mesh_factory_ptr->AddEntity(
        lf::base::RefEl::kQuad(),
        std::span<const size_type>({node[0], node[1], node[2], node[3]}),
        std::make_unique<lf::geometry::Parallelogram>(quad_coord));
  }

  // Optional: Inspect data
  // mesh_factory_ptr->PrintLists(std::cout);
  return mesh_factory_ptr->Build();
}
