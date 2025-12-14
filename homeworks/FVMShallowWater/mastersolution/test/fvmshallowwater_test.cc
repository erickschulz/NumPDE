/**
 * @file FVMShallowWater_test.cc
 * @brief NPDE homework FVMShallowWater code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../fvmshallowwater.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstddef>

/* Test in the google testing framework

  The following assertions are available, syntax
  EXPECT_XX( ....) << [anything that can be givne to std::cerr]

  EXPECT_EQ(val1, val2)
  EXPECT_NEAR(val1, val2, abs_error) -> should be used for numerical results!
  EXPECT_NE(val1, val2)
  EXPECT_TRUE(condition)
  EXPECT_FALSE(condition)
  EXPECT_GE(val1, val2)
  EXPECT_LE(val1, val2)
  EXPECT_GT(val1, val2)
  EXPECT_LT(val1, val2)
  EXPECT_STREQ(str1,str2)
  EXPECT_STRNE(str1,str2)
  EXPECT_STRCASEEQ(str1,str2)
  EXPECT_STRCASENE(str1,str2)

  "EXPECT" can be replaced with "ASSERT" when you want to program to terminate,
 if the assertion is violated.
 */

namespace FVMShallowWater::test {

TEST(FVMShallowWater, EVTest) {
  const Eigen::Vector2d u{1.0, -2.0};
  Eigen::Matrix2d DF{sweJacobian(u)};
  Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(DF);
  const Eigen::Vector2d ref_ev = eigensolver.eigenvalues().real();
  [[maybe_unused]] const Eigen::Matrix2d ref_EV =
      eigensolver.eigenvectors().real();
  auto [l1, l2] = sweLambdas(u);
  const Eigen::Vector2d lambda{l1, l2};
  EXPECT_NEAR(
      std::min((ref_ev - lambda).norm(), (ref_ev.reverse() - lambda).norm()),
      0.0, 1E-8)
      << "Mismatch of eigenvalues";
  const Eigen::Matrix2d EVs{sweEV(u)};
  Eigen::Matrix2d D = (Eigen::Matrix2d() << l1, 0.0, 0.0, l2).finished();
  EXPECT_NEAR((EVs * D * EVs.inverse() - DF).norm(), 0.0, 1E-8);
}

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix2d absMat(Eigen::Matrix2d M) {
  Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(M);
  const Eigen::Vector2d ref_ev = eigensolver.eigenvalues().cwiseAbs();
  const auto ref_EV = eigensolver.eigenvectors();
  return (ref_EV * ref_ev.asDiagonal() * ref_EV.inverse()).real();
}
/* SAM_LISTING_END_1 */

TEST(FVMShallowWater, ARabs) {
  const Eigen::Vector2d u{3.0, -2.0};
  // Compute modulus of Jacobian at state u
  auto [l1, l2] = sweLambdas(u);
  const Eigen::Matrix2d EV = sweEV(u);
  const Eigen::Matrix2d D_abs =
      (Eigen::Matrix2d() << std::abs(l1), 0.0, 0.0, std::abs(l2)).finished();
  const Eigen::Matrix2d A_abs = EV * D_abs * EV.inverse();
  // Alternative
  Eigen::Matrix2d DF{sweJacobian(u)};
  const Eigen::Matrix2d A_abs_alt{absMat(DF)};
  EXPECT_NEAR((A_abs - A_abs_alt).norm(), 0.0, 1E-8);
}

TEST(FVMShallowWater, numfluxLFSWE1) {
  Eigen::Vector2d u, v, exact_flux;
  u << 1., 2.;
  v << 1., 2.;
  exact_flux << 2., 4.5;
  Eigen::Vector2d res = numfluxLFSWE(u, v);
  ASSERT_NEAR((res - exact_flux).norm(), 0.0, 1e-8);
}

TEST(FVMShallowWater, numfluxLFSWE2) {
  Eigen::Vector2d u, v, exact_flux;
  u << 1., 2.;
  v << 3., 5.;
  exact_flux << 2., 4.5;
  Eigen::Vector2d res = numfluxLFSWE(u, v);
  ASSERT_NEAR((res - exact_flux).norm(), 0.0, 1e-8);
}

TEST(FVMShallowWater, numfluxHLLESWE1) {
  Eigen::Vector2d u, v, exact_flux;
  u << 1., 3.;
  v << .5, 4.;
  exact_flux << 3., 9.5;
  Eigen::Vector2d res = numfluxHLLESWE(u, v);
  ASSERT_NEAR((res - exact_flux).norm(), 0.0, 1e-8);
}

TEST(FVMShallowWater, numfluxHLLESWE2) {
  Eigen::Vector2d u, v, exact_flux;
  u << 1., -4.;
  v << .5, -3.;
  exact_flux << -3., 18. + 1. / 8.;
  Eigen::Vector2d res = numfluxHLLESWE(u, v);
  ASSERT_NEAR((res - exact_flux).norm(), 0.0, 1e-8);
}

TEST(FVMShallowWater, numfluxHLLESWE3) {
  Eigen::Vector2d u, v, exact_flux;
  u << 1., -4.;
  v << .5, 3.;
  exact_flux << -0.128442, 0.149371;
  Eigen::Vector2d res = numfluxHLLESWE(u, v);
  ASSERT_NEAR((res - exact_flux).norm(), 0.0, 1e-5);
}

inline double vec_norm(std::vector<double>& u) {
  const size_t N = u.size();
  double s = 0.0;
  for (int j = 0; j < N; ++j) {
    s += (u[j] * u[j]);
  }
  return std::sqrt(s);
}

TEST(FVMShallowWater, GenEvl) {
  // Simple upwind flux for linear advection
  // Data
  double a = -1.0;
  double b = 2.0;
  double T = 0.5;
  size_t N = 300;
  // Initial values
  std::vector<double> u0(N, 0.0);
  u0[100] = u0[101] = 1.0;
  // Recording intermediate approximations
  std::vector<std::vector<double>> data;
  std::vector<double> times;
  std::vector<double> u{u0};
  FVMShallowWater::FVMEvlGeneric(
      a, b, T,
      [](double h, const std::vector<double>& /*u*/) -> double { return h; }, u,
      [](double v, double /*w*/) -> double { return v; },
      [&data, &times](double t, const std::vector<double>& u) -> void {
        times.push_back(t);
        data.push_back(u);
      });
  std::cout << "FVMShallowWater::FVMEvlGeneric: " << data.size() - 1
            << " timesteps\n";
  EXPECT_NEAR(vec_norm(u0), vec_norm(u), 1E-8);
}



TEST(FVMShallowWater, RHJC1) {
  Eigen::Vector2d ul, ur;
  ul << 3., 6.;
  ur << 1., 0.36700683814454793;
  double speed;
  bool res = FVMShallowWater::checkSWERHJC(ul, ur, &speed);
  ASSERT_TRUE(res);
  EXPECT_NEAR(speed, 2.8165, 1e-3);
}

TEST(FVMShallowWater, RHJC2) {
  Eigen::Vector2d ul, ur;
  ul << 3., 6.;
  ur << 1., 1.;
  double speed;
  bool res = FVMShallowWater::checkSWERHJC(ul, ur, &speed);
  ASSERT_FALSE(res);
}

TEST(FVMShallowWater, RHJC3) {
  Eigen::Vector2d ul, ur;
  ul << 1, 0.36700683814454793;
  ur << 0.5, -0.12268279877562327;
  double speed;
  bool res = FVMShallowWater::checkSWERHJC(ul, ur, &speed);
  ASSERT_TRUE(res);
  EXPECT_NEAR(speed, 0.9793792, 1e-3);
}

TEST(FVMShallowWater, RHJC4) {
  Eigen::Vector2d ul, ur;
  ul << 1, 0.36700683814454793;
  ur << 0.5, 0.4896896369201712;
  double speed;
  bool res = FVMShallowWater::checkSWERHJC(ul, ur, &speed);
  ASSERT_TRUE(res);
  EXPECT_NEAR(speed, -0.24536559755124654, 1e-3);
}

TEST(FVMShallowWater, checkSWEPhysicalShock1) {
  Eigen::Vector2d ul, ur;
  ul << 3., 6.;
  ur << 1., 0.36700683814454793;
  bool res = FVMShallowWater::checkSWEPhysicalShock(ul, ur);
  ASSERT_TRUE(res);
}

TEST(FVMShallowWater, checkSWEPhysicalShock2) {
  Eigen::Vector2d ul, ur;
  ul << 3., 6.;
  ur << 1., 3.632993161855452;
  bool res = FVMShallowWater::checkSWEPhysicalShock(ul, ur);
  ASSERT_FALSE(res);
}

TEST(FVMShallowWater, checkSWEPhysicalShock3) {
  Eigen::Vector2d ul, ur;
  ul << 1, 0.36700683814454793;
  ur << 0.5, -0.12268279877562327;
  double speed;
  bool res = FVMShallowWater::checkSWEPhysicalShock(ul, ur);
  ASSERT_TRUE(res);
}

TEST(FVMShallowWater, checkSWEPhysicalShock4) {
  Eigen::Vector2d ul, ur;
  ul << 1, 0.36700683814454793;
  ur << 0.5, 0.4896896369201712;
  bool res = FVMShallowWater::checkSWEPhysicalShock(ul, ur);
  ASSERT_FALSE(res);
}


TEST(FVMShallowWater, B_isPhysicalTwoShockSolution1) {
  Eigen::Vector2d ul, us, ur;
  ul << 0.5, 0.4896896369201712;
  us << 1, 0.36700683814454793;
  ur << 0.5, -0.12268279877562327;

  bool res = FVMShallowWater::isPhysicalTwoShockSolution(ul, us, ur);
  ASSERT_TRUE(res);
}

TEST(FVMShallowWater, B_isPhysicalTwoShockSolution2) {
  Eigen::Vector2d ul, us, ur;
  ul << 0.5, 0.4896896369201712;
  us << 1, 0.36700683814454793;
  ur << 0.5, 0.4896896369201712;

  bool res = FVMShallowWater::isPhysicalTwoShockSolution(ul, us, ur);
  ASSERT_FALSE(res);
}

TEST(FVMShallowWater, solveSWE) {
  // Simple upwind flux for linear advection
  // Data
  double a = -1.0;
  double b = 2.0;
  double T = 0.5;
  size_t N = 300;
  // Initial values
  auto u0 = [](double x) -> Eigen::Vector2d {
    if (x > 0 && x < 1)
      return (Eigen::Vector2d() << 1. - x * (1. - x), 1.).finished();

    return (Eigen::Vector2d() << 1., 0.).finished();
  };
  // Recording intermediate approximations
  auto res =
      FVMShallowWater::solveSWE(T, N, &FVMShallowWater::numfluxLFSWE, u0);
  Eigen::Vector2d exact_at_200;
  exact_at_200 << 0.411047, 0.294142;
  ASSERT_NEAR((res[200] - exact_at_200).norm(), 0., 1e-4);
}

}  // namespace FVMShallowWater::test
