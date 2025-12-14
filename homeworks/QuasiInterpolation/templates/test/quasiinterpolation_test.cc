/**
 * @file quasiinterpolation_test.cc
 * @brief NPDE exam QuasiInterpolation code
 * @author Oliver Rietmann
 * @date 15.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../quasiinterpolation.h"

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <memory>
#include <utility>

namespace QuasiInterpolation::test {

TEST(QuasiInterpolation, findKp) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  lf::mesh::utils::CodimMeshDataSet<
      std::pair<const lf::mesh::Entity*, unsigned int>>
      KpMeshDataSet = findKp(mesh_p);

  unsigned int vertexIndex = 9;       // 8
  unsigned int localVertexIndex = 1;  // 0
  unsigned int triangleIndex = 12;    // 11

  const lf::mesh::Entity* vertex = mesh_p->EntityByIndex(2, vertexIndex);
  std::pair<const lf::mesh::Entity*, unsigned int> Kp = KpMeshDataSet(*vertex);

  EXPECT_EQ(localVertexIndex, Kp.second);
  EXPECT_NE(Kp.first, nullptr);
  if (Kp.first != nullptr) {
    EXPECT_EQ(triangleIndex, mesh_p->Index(*Kp.first));
  }
}

TEST(QuasiInterpolation, quasiInterpolate) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  lf::uscalfe::FeSpaceLagrangeO1<double> fe_space(mesh_p);

  // For polynomials of degree 1, the quasi projection yields the exact result
  auto f = [](Eigen::Vector2d x) -> double {
    return 1.5 * x(0) - 2.0 * x(1) + 0.5;
  };
  lf::mesh::utils::MeshFunctionGlobal mf(f);

  Eigen::VectorXd coefficients =
      QuasiInterpolation::quasiInterpolate(fe_space, mf);
  Eigen::VectorXd coefficients_ref = lf::fe::NodalProjection(fe_space, mf);

  double tol = 1.0e-12;
  double error = (coefficients - coefficients_ref).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
}

}  // namespace QuasiInterpolation::test
