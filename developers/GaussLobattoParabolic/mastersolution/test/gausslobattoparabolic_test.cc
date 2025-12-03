/**
 * @file gausslobattoparabolic_test.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../gausslobattoparabolic.h"

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <memory>

namespace GaussLobattoParabolic::test {

Eigen::Matrix<double, 10, 10> getM() {
  int N = 10;
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);

  M.row(0) << 0.888888888888889, 0.0902777777777778, 0.131944444444444,
      0.0277777777777778, 0.0, 0.0, 0.0416666666666667, 0.215277777777778, 0.25,
      0.1875;
  M.row(1) << 0.0902777777777778, 0.375, 0.0833333333333334, 0.111111111111111,
      0.104166666666667, 0.0, 0.0, 0.0, 0.0, 0.0277777777777778;
  M.row(2) << 0.131944444444444, 0.0833333333333334, 0.625, 0.0,
      0.104166666666667, 0.145833333333333, 0.159722222222222,
      0.0416666666666667, 0.0, 0.0;

  return M;
}

Eigen::Matrix<double, 10, 10> getA() {
  int N = 10;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);

  A.row(0) << 4.16763791763792, -1.11363636363636, -0.637820512820513,
      0.0757575757575758, 0.0, 0.0, -0.64957264957265, -0.637820512820513,
      -1.08333333333333, -0.121212121212121;
  A.row(1) << -1.11363636363636, 4.29545454545455, -0.75, -0.363636363636364,
      -0.75, 0.0, 0.0, 0.0, 0.0, -1.31818181818182;
  A.row(2) << -0.637820512820513, -0.75, 3.58173076923077, 0.0, -0.75, -0.75,
      -0.525641025641026, -0.168269230769231, 0.0, 0.0;

  for (int i = 3; i < N; ++i) A(i, i) = 1.0;

  return A;
}

TEST(GaussLobattoParabolic, initMbig) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::MatrixXd M = initMbig(fe_space).makeDense();

  int N = 10;
  EXPECT_EQ(N, M.rows());
  EXPECT_EQ(N, M.cols());

  if (N == M.rows() && N == M.cols()) {
    Eigen::MatrixXd M_ref = getM();

    double tol = 1.0e-8;
    double error = (M - M_ref).lpNorm<Eigen::Infinity>();
    ASSERT_NEAR(0.0, error, tol);
  }
}

TEST(GaussLobattoParabolic, initAbig) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::MatrixXd A = initAbig(fe_space).makeDense();

  int N = 10;
  EXPECT_EQ(N, A.rows());
  EXPECT_EQ(N, A.cols());

  if (N == A.rows() && N == A.cols()) {
    Eigen::MatrixXd A_ref = getA();

    double tol = 1.0e-8;
    double error = (A - A_ref).lpNorm<Eigen::Infinity>();
    ASSERT_NEAR(0.0, error, tol);
  }
}

TEST(GaussLobattoParabolic, RHSProvider) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  RHSProvider rhs_provider(fe_space->LocGlobMap(),
                           (double (*)(double)) & std::sin);

  Eigen::VectorXd phi = rhs_provider(1.0);
  int N = phi.size();

  int N_ref = 10;
  EXPECT_EQ(N, N_ref);

  if (N == N_ref) {
    Eigen::VectorXd phi_ref(N_ref);
    phi_ref << 0.0, 0.0, 0.0, 0.841470984807897, 0.841470984807897,
        0.841470984807897, 0.841470984807897, 0.841470984807897,
        0.841470984807897, 0.841470984807897;

    double tol = 1.0e-8;
    double error = (phi - phi_ref).lpNorm<Eigen::Infinity>();
    ASSERT_NEAR(0.0, error, tol);
  }
}

TEST(GaussLobattoParabolic, evolveIBVPGaussLobatto) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  double T = 1.0;
  unsigned int M = 100;
  auto g = [T](double t) {
    const double PI = 3.14159265358979323846;
    return t < 1.0 ? std::sin(0.5 * PI * t) : 1.0;
  };
  Eigen::VectorXd mu = evolveIBVPGaussLobatto(fe_space, T, M, g);

  int N = 10;
  Eigen::VectorXd mu_ref(N);
  mu_ref << 0.606154338132087, 0.732422646475846, 0.655377750264433, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0;

  double tol = 1.0e-6;
  double error = (mu - mu_ref).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
}

}  // namespace GaussLobattoParabolic::test
