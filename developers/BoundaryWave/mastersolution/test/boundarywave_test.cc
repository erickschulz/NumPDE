/**
 * @file boundarywave_test.cc
 * @brief NPDE homework "BoundaryWave" code
 * @author Philipp Lindenberger
 * @date 03.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <memory>
#include <utility>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "../boundarywave.h"

namespace BoundaryWave::test {

auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/simple.msh");
auto mesh_p = reader.mesh();

TEST(BoundaryWave, buildM) {
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  lf::assemble::COOMatrix M_coo = buildM(fe_space_p);

  Eigen::MatrixXd M = M_coo.makeDense();

  // std::cout << M;
  double eps = 1.0e-5;
  Eigen::Matrix<double, 5, 5> reference_M;

  reference_M << 0.666667, 0.166667, 0, 0.166667, 0, 0.166667, 0.666667,
      0.166667, 0, 0, 0, 0.166667, 0.666667, 0.166667, 0, 0.166667, 0, 0.166667,
      0.666667, 0, 0, 0, 0, 0, 0;
  ASSERT_EQ(reference_M.size(), M.size());
  ASSERT_NEAR((reference_M - M).lpNorm<Eigen::Infinity>(), 0.0, eps);
}

TEST(BoundaryWave, buildA) {
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  lf::assemble::COOMatrix A_coo = buildA(fe_space_p);
  Eigen::MatrixXd A = A_coo.makeDense();

  double eps = 1.0e-5;
  Eigen::Matrix<double, 5, 5> reference_A;

  reference_A << 1.33333, 0, 0, 0, -1.33333, 0, 1.66667, 0, 0, -1.66667, 0, 0,
      2, 0, -2, 0, 0, 0, 1.66667, -1.66667, -1.33333, -1.66667, -2, -1.66667,
      6.66667;
  ASSERT_EQ(reference_A.size(), A.size());
  ASSERT_NEAR((reference_A - A).lpNorm<Eigen::Infinity>(), 0.0, eps);
}

TEST(BoundaryWave, InterpolateInitialData) {
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0 = [](const Eigen::Vector2d& x) -> double { return x[0]; };
  auto v0 = [](const Eigen::Vector2d& x) -> double { return x[1]; };

  std::pair<Eigen::VectorXd, Eigen::VectorXd> initialData =
      interpolateInitialData(fe_space_p, std::move(u0), std::move(v0));

  Eigen::VectorXd reference_u0(5);
  reference_u0 << 0, 1, 1, 0, 0.5;

  Eigen::VectorXd reference_v0(5);
  reference_v0 << 0, 0, 1, 1, 0.5;

  double eps = 1.0e-5;
  ASSERT_TRUE(initialData.first.size() == reference_u0.size());
  ASSERT_TRUE(initialData.second.size() == reference_v0.size());
  for (int i = 0; i < 5; i++) {
    ASSERT_NEAR(reference_u0(i), initialData.first(i), eps);
    ASSERT_NEAR(reference_v0(i), initialData.second(i), eps);
  }
}

TEST(BoundaryWave, solveBoundaryWave) {
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0 = [](const Eigen::Vector2d& x) -> double { return x[0]; };
  auto v0 = [](const Eigen::Vector2d& x) -> double { return x[1]; };

  Eigen::VectorXd discrete_solution =
      solveBoundaryWave(fe_space_p, u0, v0, 1.0, 100);

  Eigen::VectorXd reference_solution(5);

  reference_solution << 0.584002, 0.767904, 1.23746, 1.41064, 0.982673;

  double eps = 1.0e-5;
  ASSERT_EQ(reference_solution.rows(), discrete_solution.rows());
  ASSERT_NEAR(
      (reference_solution - discrete_solution).lpNorm<Eigen::Infinity>(), 0.0,
      eps);
}
}  // namespace BoundaryWave::test
