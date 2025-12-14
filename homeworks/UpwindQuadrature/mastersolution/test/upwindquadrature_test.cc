/**
 * @file upwindquadrature_test.cc
 * @brief NPDE homework UpwindQuadrature code
 * @author  Philippe Peter
 * @date June 2020
 * @copyright Developed at ETH Zurich
 */

#include "../upwindquadrature.h"

#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <memory>

#include "../../../../lecturecodes/ConvectionDiffusion/convection_emp.h"

namespace UpwindQuadrature::test {

TEST(UpwindQuadrature, upwind_convection_element_matrix_provider_1) {
  // construct a triangular tensor product mesh on the unit square
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(2)
      .setNumYCells(2);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

  // initialize masses and velocity field
  auto masses = UpwindQuadrature::initializeMasses(mesh_p);
  const auto v = [](const Eigen::Vector2d& /*x*/) {
    return (Eigen::Vector2d() << 1, 2).finished();
  };

  // initialize already implemented reference provider and
  // the upwind provider
  ConvectionDiffusion::ConvectionElementMatrixProvider reference(v);
  UpwindConvectionElementMatrixProvider upwind(v, masses);

  const lf::mesh::Entity& element_0 = *(mesh_p->EntityByIndex(0, 0));
  const lf::mesh::Entity& element_1 = *(mesh_p->EntityByIndex(0, 1));

  // element 0 is at nonoe of the corners the upwind triangle.
  EXPECT_NEAR(upwind.Eval(element_0).norm(), 0.0, 1E-15);

  // element 1 is the upwind triangle at corner 2.
  Eigen::MatrixXd reference_eval = reference.Eval(element_1);
  Eigen::MatrixXd upwind_eval = upwind.Eval(element_1);

  EXPECT_NEAR(upwind_eval.row(0).norm(), 0.0, 1E-14);
  EXPECT_NEAR(upwind_eval.row(1).norm(), 0.0, 1E-14);
  EXPECT_NEAR((upwind_eval.row(2) - 6.0 * reference_eval.row(2)).norm(), 0.0,
              1E-14);
}

TEST(UpwindQuadrature, upwind_convection_element_matrix_provider_2) {
  // construct a triangular tensor product mesh on the unit square
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(2)
      .setNumYCells(2);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

  // initialize masses and velocity field
  auto masses = UpwindQuadrature::initializeMasses(mesh_p);
  const auto v = [](const Eigen::Vector2d& /*x*/) {
    return (Eigen::Vector2d() << -2, -1).finished();
  };

  // initialize already implemented reference provider and
  // the upwind provider
  ConvectionDiffusion::ConvectionElementMatrixProvider reference(v);
  UpwindConvectionElementMatrixProvider upwind(v, masses);

  const lf::mesh::Entity& element_1 = *(mesh_p->EntityByIndex(0, 1));

  // element 1 is the upwind triangle at corner 0.
  Eigen::MatrixXd reference_eval = reference.Eval(element_1);
  Eigen::MatrixXd upwind_eval = upwind.Eval(element_1);

  EXPECT_NEAR((upwind_eval.row(0) - 2.0 * reference_eval.row(0)).norm(), 0.0,
              1E-14);
  EXPECT_NEAR(upwind_eval.row(1).norm(), 0.0, 1E-14);
  EXPECT_NEAR(upwind_eval.row(2).norm(), 0.0, 1E-14);
}

TEST(UpwindQuadrature, upwind_convection_element_matrix_provider_3) {
  // construct a triangular tensor product mesh on the unit square
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(2)
      .setNumYCells(2);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

  // initialize masses and velocity field
  auto masses = UpwindQuadrature::initializeMasses(mesh_p);
  const auto v = [](const Eigen::Vector2d& /*x*/) {
    return (Eigen::Vector2d() << 0, -1).finished();
  };

  // initialize already implemented reference provider and
  // the upwind provider
  ConvectionDiffusion::ConvectionElementMatrixProvider reference(v);
  UpwindConvectionElementMatrixProvider upwind(v, masses);

  const lf::mesh::Entity& element_1 = *(mesh_p->EntityByIndex(0, 1));

  // at corner 1 of element 1, -v(a^1) points along the edge
  //--> contribution split between triangle sharing that edge.
  Eigen::MatrixXd reference_eval = reference.Eval(element_1);
  Eigen::MatrixXd upwind_eval = upwind.Eval(element_1);

  EXPECT_NEAR(upwind_eval.row(0).norm(), 0.0, 1E-14);
  EXPECT_NEAR((upwind_eval.row(1) - 1.5 * reference_eval.row(1)).norm(), 0.0,
              1E-14);
  EXPECT_NEAR(upwind_eval.row(2).norm(), 0.0, 1E-14);
}

}  // namespace UpwindQuadrature::test
