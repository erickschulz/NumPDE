/**
 * @file neumanndatarecovery_test.cc
 * @brief NPDE homework NeumannDataRecovery code
 * @author R. Hiptmair
 * @date July 7, 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../neumanndatarecovery.h"

#include <gtest/gtest.h>

#include "lf/mesh/test_utils/test_meshes.h"

namespace NeumannDataRecovery::test {
TEST(NeumannDataRecovery, UnitNormals) {
  // Generate a test mesh
  const int selector = 0;  // Hybrid mesh
  const double scale = 1.0 / 3.0;
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector, scale);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Loop through all cells of the mesh
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Obtain geometry information
    const lf::geometry::Geometry* geo_ptr = cell->Geometry();
    // Compute exterior unit normals
    const Eigen::Matrix<double, 2, 4> unit_normals =
        exteriorUnitNormals(*geo_ptr);
    // Run through the edges of the current cell and check whether the normals
    // are really orthogonal to the edge direction vectors
    std::span<const lf::mesh::Entity* const> edges{cell->SubEntities(1)};
    int ed_cnt = 0;
    for (const lf::mesh::Entity* edge : edges) {
      LF_ASSERT_MSG((edge->RefEl() == lf::base::RefEl::kSegment()),
                    "Edge must be a segment");
      const Eigen::Matrix<double, 2, 2> ed_vt =
          lf::geometry::Corners(*(edge->Geometry()));
      // std::cout << "Edge " << ed_cnt << " = \n" << ed_vt << std::endl;
      const Eigen::Vector2d edv = ed_vt.col(1) - ed_vt.col(0);
      EXPECT_NEAR(unit_normals.col(ed_cnt).norm(), 1.0, 1.0E-10);
      EXPECT_NEAR(edv.dot(unit_normals.col(ed_cnt)), 0.0, 1E-10);
      ed_cnt++;
    }
  }
}

TEST(NeumannDataRecovery, PwconstNeumannData) {
  // Generate a test mesh of the unit square
  const int selector = 3;  // Mesh comprising only triangles
  const double scale = 1.0 / 3.0;
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector, scale);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Interpolate a linear function into the finite element space
  // This function can be represented exactly.
  lf::mesh::utils::MeshFunctionGlobal mf_f{[](Eigen::Vector2d x) -> double {
    return 1.0 + 3.0 * x[0] - 2.0 * x[1];
  }};
  const Eigen::VectorXd f_dof = lf::fe::NodalProjection(*fe_space, mf_f);
  lf::mesh::utils::CodimMeshDataSet<double> neu_dat =
      getNeumannData(fe_space, f_dof);
  // Compute integral of Neumann data. To that end run through all edges on the
  // boundary and add edge contributions Obtain flags indicating edges on the
  // boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bded_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  // Loop over all edges of the mesh
  double s = 0.0;
  bool some_nonzero = false;
  for (const lf::mesh::Entity* edge : fe_space->Mesh()->Entities(1)) {
    LF_ASSERT_MSG(edge->RefEl() == lf::base::RefEl::kSegment(),
                  "Edge must be of segment type!");
    if (bded_flags(*edge)) {
      const double neudat_ptval = neu_dat(*edge);
      s += lf::geometry::Volume(*(edge->Geometry())) * neudat_ptval;
      if (std::fabs(neudat_ptval) > 1e-8) {
        some_nonzero = true;
      }
    }
  }
  EXPECT_NEAR(s, 0.0, 1E-6);
  ASSERT_TRUE(some_nonzero);
}
}  // namespace NeumannDataRecovery::test
