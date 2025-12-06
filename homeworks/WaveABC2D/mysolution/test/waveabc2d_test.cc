/**
 * @file waveabc2d_test.cc
 * @brief NPDE homework "WaveABC2D" code
 * @author Philipp Lindenberger
 * @date 03.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

// Lehrfem++ includes
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include "../waveabc2d.h"

namespace WaveABC2D::test {

TEST(WaveABC2D, scalarImplicitTimestepping) {
  double eps = 1.0e-5;
  double epsilon = 0.5;
  int M = 10;
  Eigen::VectorXd reference_solution(M + 1);
  reference_solution << 1.0, 1.08531, 1.14196, 1.17159, 1.1762, 1.15807,
      1.11971, 1.06376, 0.992949, 0.910024, 0.817708;

  Eigen::VectorXd student_solution = scalarImplicitTimestepping(epsilon, M);

  ASSERT_NEAR((reference_solution - student_solution).lpNorm<Eigen::Infinity>(),
              0.0, eps);
}

// Test with constant initial data
TEST(WaveABC2D, WaveABC2DTimestepper_const) {
  double eps = 1.0e-5;
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto rho = [](Eigen::Vector2d) -> double { return 4.0; };
  auto mu0 = [](const Eigen::Vector2d& x) -> double { return 1.0; };
  auto nu0 = [](const Eigen::Vector2d& x) -> double { return 1.0; };
  auto stepper =
      WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)>(
          fe_space_p, rho, 500, 1.0);
  auto student_solution = stepper.solveWaveABC2D(mu0, nu0);

  Eigen::VectorXd reference_solution(13);
  reference_solution << 1.57776, 1.70942, 1.57249, 2.11247, 2.1447, 1.67896,
      1.71236, 2.14945, 2.11761, 1.64201, 1.57128, 1.74306, 1.62853;

  ASSERT_EQ(student_solution.size(), reference_solution.size());
  ASSERT_NEAR((reference_solution - student_solution).lpNorm<Eigen::Infinity>(),
              0.0, eps);
}

// Test with original initial data
TEST(WaveABC2D, WaveABC2DTimestepper) {
  double eps = 1.0e-5;
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto rho = [](Eigen::Vector2d) -> double { return 4.0; };

  auto mu0 = [](const Eigen::Vector2d& x) -> double {
    return std::sin(x.norm());
  };
  auto nu0 = [](const Eigen::Vector2d& x) -> double { return std::cos(x(1)); };

  auto stepper =
      WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)>(
          fe_space_p, rho, 500, 1.0);

  Eigen::VectorXd student_solution = stepper.solveWaveABC2D(mu0, nu0);

  Eigen::VectorXd reference_solution(13);
  reference_solution << 0.97316, 1.50735, 0.875053, 1.47574, 1.17419, 0.613036,
      0.883386, 0.342844, 0.0600508, -0.343508, -0.232577, -0.621964, -1.3619;

  ASSERT_EQ(student_solution.size(), reference_solution.size());
  ASSERT_NEAR((reference_solution - student_solution).lpNorm<Eigen::Infinity>(),
              0.0, eps);
}

TEST(WaveABC2D, energies) {
  double eps = 1.0e-4;
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto rho = [](Eigen::Vector2d) -> double { return 4.0; };

  auto mu0 = [](const Eigen::Vector2d& x) -> double {
    return std::sin(x.norm());
  };
  auto nu0 = [](const Eigen::Vector2d& x) -> double { return std::cos(x(1)); };

  auto stepper =
      WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)>(
          fe_space_p, rho, 500, 1.0);

  Eigen::VectorXd sol = stepper.solveWaveABC2D(mu0, nu0);
  double student_solution = stepper.energies();

  double reference_solution = 11.4534;

  ASSERT_NEAR(reference_solution - student_solution, 0, eps);
}

}  // namespace WaveABC2D::test
