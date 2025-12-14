/**
 * @file leapfrogdissipativewave_test.cc
 * @brief NPDE homework LeapfrogDissipativeWave code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../leapfrogdissipativewave.h"

#include <gtest/gtest.h>

#include <iostream>

#include "lf/mesh/test_utils/test_meshes.h"

namespace LeapfrogDissipativeWave::test {

TEST(LeapfrogDissipativeWave, timestepDissipativeWaveEquation) {
  // Generate a triangular test mesh of the unit square
  const int selector = 3;  // Mesh comprising only triangles
  const double scale = 1.0 / 3.0;
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector, scale);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Set up global FE space; quadratic Lagrangian finite elements
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Initial Conditions
  lf::mesh::utils::MeshFunctionGlobal mf_u0{[](Eigen::Vector2d x) -> double {
    return std::cos(M_PI * (2 * x[0] - 1)) * std::cos(M_PI * (2 * x[1] - 1));
  }};
  lf::mesh::utils::MeshFunctionConstant mf_v0{1.0};
  // Interpolate initial data into FE space
  const Eigen::VectorXd u0_vec = lf::fe::NodalProjection(*fe_space_p, mf_u0);
  const Eigen::VectorXd v0_vec = lf::fe::NodalProjection(*fe_space_p, mf_v0);
  // Propagate initial conditions in time
  const Eigen::VectorXd sol =
      timestepDissipativeWaveEquation(fe_space_p, 1, 100, u0_vec, v0_vec);
  // Correct solution
  Eigen::VectorXd ref_sol(40);
  ref_sol << 0.157742, 0.229119, 0.195882, 0.329603, 0.289977, 0.210851,
      0.178861, 0.317668, 0.27685, 0.220541, 0.161492, 0.20694, 0.179893,
      0.20649, 0.256688, 0.185787, 0.24507, 0.312891, 0.27881, 0.253802,
      0.209621, 0.323596, 0.318452, 0.271053, 0.305812, 0.285747, 0.284253,
      0.304414, 0.215171, 0.313654, 0.225985, 0.278874, 0.284611, 0.310179,
      0.279783, 0.259143, 0.293064, 0.217791, 0.232835, 0.225115;
  ASSERT_NEAR((sol - ref_sol).lpNorm<Eigen::Infinity>(), 0, 1e-6);
}

}  // namespace LeapfrogDissipativeWave::test
