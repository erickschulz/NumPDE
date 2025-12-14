/**
 * @file
 * @brief Test for solving source-free heat equation and computation of H1
 *  	  seminorms on different triangular meshes and refinement levels
 * @author Julien Gacon
 * @date   March 2019
 * @copyright MIT License
 */

#include "../unstablebvp.h"

// General includes
#include <cmath>
#include <memory>
// Google Test
#include <gtest/gtest.h>
// Lehrfempp
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

TEST(UnstableBVP, TopMesh) {
  // Define the number of refinement levels we want for our mesh
  const int reflevels = 7;

  // Get a hierachy of refined meshes
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      UnstableBVP::createMeshHierarchy(reflevels, "top");
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  for (int level = 0; level <= reflevels; ++level) {
    // Get the mesh pointer
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);

    // Get the seminorm
    double h1 = UnstableBVP::solveTemperatureDistribution(mesh_p);

    EXPECT_NEAR(h1, 0., 1e-15);
  }
}

TEST(UnstableBVP, BottomMesh) {
  // Define the number of refinement levels we want for our mesh
  const int reflevels = 7;

  // Get a hierachy of refined meshes
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      UnstableBVP::createMeshHierarchy(reflevels, "bottom");
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  for (int level = 0; level <= reflevels; ++level) {
    // Get the mesh pointer
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);

    // Get the seminorm
    double h1 = UnstableBVP::solveTemperatureDistribution(mesh_p);

    EXPECT_NEAR(h1, 0.707107, 1e-5);
  }
}

TEST(UnstableBVP, CenterMesh) {
  // Get a hierachy of refined meshes with 7 reflevels
  // levels are hardcoded for they shouldn't be changed due to dependencies
  // in the test
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      UnstableBVP::createMeshHierarchy(7, "center");
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  // Get the seminorm of level 7 mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(7);
  double h1_uL = UnstableBVP::solveTemperatureDistribution(mesh_p);

  for (int level : {0, 3, 6}) {
    // Get the mesh pointer
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);

    // Get the seminorm
    double h1 = UnstableBVP::solveTemperatureDistribution(mesh_p);

    switch (level) {
      case 0:
        EXPECT_NEAR(h1, 1.06066, 1e-5);
        EXPECT_NEAR(std::abs(h1 - h1_uL), 0.967539, 1e-6);
        break;
      case 3:
        EXPECT_NEAR(h1, 1.50899, 1e-5);
        EXPECT_NEAR(std::abs(h1 - h1_uL), 0.519206, 1e-6);
        break;
      case 6:
        EXPECT_NEAR(h1, 1.913105, 1e-6);
        EXPECT_NEAR(std::abs(h1 - h1_uL), 0.115094, 1e-6);
        break;
      default:
        throw "Test for this level not implemented.";
    }
  }
}
