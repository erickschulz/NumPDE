/**
 * @file extendedmuscl_test_mastersolution.cc
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include "../extendedmuscl.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>

#include "../slopelimfluxdiff.h"

namespace ExtendedMUSCL::test {

constexpr double PI = 3.14159265358979323846;

TEST(ExtendedMUSCL, logGodunovFlux) {
  // setting
  Eigen::VectorXd v(13);
  v << 4.22, 3.82, 2.39, 2.84, 4.60, 1.88, 4.23, 3.14, 0.37, 0.3, 0.1, 0.5, 1.5;
  Eigen::VectorXd w(13);
  w << 4.63, 3.89, 1.81, 4.80, 2.05, 1.81, 0.63, 4.08, 0.40, 0.1, 0.3, 0.5, 1.5;

  // my solution
  Eigen::VectorXd F_GD = v.binaryExpr(w, std::ref(logGodunovFlux));

  // reference
  Eigen::VectorXd F_GD_ref(13);
  F_GD_ref << 1.85610424036222, 1.29975661440261, -0.307608855395228,
      0.124403508171646, 2.41985899607723, -0.693209059537307, 1.87051443063613,
      0.452859591749309, -0.766516292749662, -0.330258509299405,
      -0.661191841297781, -0.846573590279973, -0.891802337837753;

  // compare
  double tol = 1.0e-8;
  double error = (F_GD_ref - F_GD).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(0.0, error, tol);
}

TEST(ExtendedMUSCL, limiterMC) {
  // setting
  int N = 10;
  Eigen::VectorXd mu_left(N);
  mu_left << 4.22, 3.82, 2.39, 2.84, 4.60, 1.88, 4.23, 3.14, 0.37, 4.64;
  Eigen::VectorXd mu_center(N);
  mu_center << 4.63, 3.89, 1.81, 4.80, 2.05, 1.81, 0.63, 4.08, 0.40, 0.99;
  Eigen::VectorXd mu_right(N);
  mu_right << 0.33, 3.22, 1.80, 0.57, 1.46, 0.92, 2.32, 0.26, 1.30, 3.05;

  // my solution
  Eigen::VectorXd scaled_slope(N);
  for (int n = 0; n < N; ++n)
    scaled_slope(n) = limiterMC(mu_left(n), mu_center(n), mu_right(n));

  // reference
  Eigen::VectorXd scaled_slope_ref(N);
  scaled_slope_ref << 0.0, 0.0, -0.02, 0.0, -1.18, -0.14, 0.0, 0.0, 0.06, 0.0;

  // compare
  double tol = 1.0e-8;
  double error = (scaled_slope_ref - scaled_slope).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(0.0, error, tol);
}

TEST(ExtendedMUSCL, slopelimfluxdiffper) {
  // setting
  unsigned int n = 20;
  auto central_slope = [](double mu_left, double mu_center, double mu_right) {
    return 0.5 * (mu_right - mu_left);
  };
  auto f = [](double u) { return 0.5 * u * u; };  // flux from Burger's equation
  auto central_flux = [&f](double v, double w) { return 0.5 * (f(v) + f(w)); };
  auto u = [](double x) { return std::sin(2.0 * PI * x) + 2.0; };
  double h = 1.0 / n;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0.5 * h, 1.0 - 0.5 * h);
  Eigen::VectorXd mu = x.unaryExpr(u);

  // my solution
  Eigen::VectorXd rhs = slopelimfluxdiffper(mu, central_flux, central_slope);

  // reference
  Eigen::VectorXd rhs_ref(n);
  rhs_ref << 0.674247307276668, 0.692128672550765, 0.605903098384332,
      0.415428324537427, 0.147931969198225, -0.147931969198225,
      -0.415428324537427, -0.605903098384332, -0.692128672550764,
      -0.674247307276668, -0.576478946599618, -0.436167781273887,
      -0.289518037184678, -0.159467433260549, -0.0501636085211743,
      0.0501636085211742, 0.159467433260548, 0.289518037184678,
      0.436167781273886, 0.576478946599618;

  // compare
  double tol = 1.0e-6;
  double error = (rhs_ref - rhs).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(0.0, error, tol);
}

TEST(ExtendedMUSCL, sspEvolop) {
  // setting
  Eigen::MatrixXd A(3, 3);
  A << 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
  auto f = [A](Eigen::Vector3d y) { return A * y; };
  const Eigen::Vector3d y0 = {0.333, -0.75, 0.9};
  double T = PI;
  int n = 200;
  double tau = T / n;

  // my solution
  Eigen::Vector3d y = y0;
  for (int i = 0; i < n; ++i) y = sspEvolop(f, y, tau);

  // reference
  Eigen::MatrixXd expTA(3, 3);
  expTA << std::cos(T), -std::sin(T), 0.0, std::sin(T), std::cos(T), 0.0, 0.0,
      0.0, std::exp(-T);
  Eigen::Vector3d y_ref = expTA * y0;

  // compare
  double error = (y_ref - y).squaredNorm();
  double tol = 1.0e-6;
  EXPECT_NEAR(0.0, error, tol);
}

TEST(ExtendedMUSCL, solveClaw) {
  // setting
  auto u0 = [](double x) { return 0.25 < x && x < 0.75 ? 2.0 : 1.0; };
  double T = 1.0;
  unsigned int n = 20;

  // my solution
  Eigen::VectorXd muT = solveClaw(u0, T, n);

  // reference
  Eigen::VectorXd muT_ref(n);
  muT_ref << 1.99885475867196, 2.00304160649059, 1.74443725020065,
      1.06007730386849, 1, 1.04277512969966, 1.09389708054387, 1.14865408584807,
      1.20623124804183, 1.26648053287162, 1.32935420687853, 1.39479385511108,
      1.4627046662712, 1.53292720003324, 1.60519872043032, 1.67900386122176,
      1.75335347818965, 1.82738742724492, 1.8973284584603, 1.95349912992223;

  // compare
  double tol = 1.0e-2;
  double error = (muT_ref - muT).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(0.0, error, tol);
}

}  // namespace ExtendedMUSCL::test
