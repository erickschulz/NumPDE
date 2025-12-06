/**
 * @file clempiricflux_test_mastersolution.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>

// In the interest of not changing the problem text over at
// https://gitlab.math.ethz.ch/ralfh/npdeflipped, I will tell clang-tidy to
// ignore the bugprone-suspicious-include warning. (Manuel Saladin, 2024-05-28)
#include "../clempiricflux.h"
#include "../solvecauchyproblem.h"
#include "../uniformcubicspline.h"

namespace CLEmpiricFlux::test {

constexpr double PI = 3.14159265358979323846;

TEST(CLEmpiricFlux, UniformCubicSpline_getJ) {
  double a = -1.0;
  double b = 2.0;
  unsigned int n = 3;

  Eigen::VectorXd x(10);
  Eigen::VectorXi j(10);
  x << -1.0, -0.9, -0.1, 0.0, 0.1, 0.9, 1.0, 1.1, 1.9, 2.0;
  j << 1, 1, 1, 2, 2, 2, 3, 3, 3, 3;

  double tol = 1.e-12;
  for (int i = 0; i < x.size(); ++i) {
    int j_test = getJ(a, b, n, x(i));
    EXPECT_NEAR(j(i), j_test, tol);
  }
}

TEST(CLEmpiricFlux, UniformCubicSpline_class) {
  int a = -1.0;
  int b = 2.0;
  unsigned int n = 5;

  Eigen::Vector4d c = {0.333, -0.5, 1.5, -0.666};
  auto f_lambda = [&c](double x) {
    return c(0) + c(1) * x + c(2) * Square(x) + c(3) * Cube(x);
  };
  auto df_lambda = [&c](double x) {
    return c(1) + 2.0 * c(2) * x + 3.0 * c(3) * Square(x);
  };
  auto M_lambda = [&c](double x) { return 2.0 * c(2) + 6.0 * c(3) * x; };

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n + 1, a, b);
  Eigen::VectorXd f = x.unaryExpr(f_lambda);
  Eigen::VectorXd M = x.unaryExpr(M_lambda);
  UniformCubicSpline spline(a, b, f, M);

  unsigned int m = 50;
  Eigen::VectorXd x_ref = Eigen::VectorXd::LinSpaced(m, a, b);
  double tol = 1.e-12;
  double error;

  // Test operator()
  Eigen::VectorXd f_ref = x_ref.unaryExpr(f_lambda);
  Eigen::VectorXd f_approx = x_ref.unaryExpr(spline);
  error = (f_ref - f_approx).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(error, 0.0, tol);

  // Test derivative
  Eigen::VectorXd df_ref = x_ref.unaryExpr(df_lambda);
  auto wrapped = [&spline](double u) { return spline.derivative(u); };
  Eigen::VectorXd df_approx = x_ref.unaryExpr(wrapped);
  error = (df_ref - df_approx).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(CLEmpiricFlux, findSupport) {
  int a = -1.0;
  int b = 2.0;
  unsigned int n = 5;

  Eigen::Vector4d c = {0.333, -0.5, 1.5, -0.666};
  auto f_lambda = [&c](double x) {
    return c(0) + c(1) * x + c(2) * Square(x) + c(3) * Cube(x);
  };
  auto M_lambda = [&c](double x) { return 2.0 * c(2) + 6.0 * c(3) * x; };

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n + 1, a, b);
  Eigen::VectorXd f = x.unaryExpr(f_lambda);
  Eigen::VectorXd M = x.unaryExpr(M_lambda);
  UniformCubicSpline spline(a, b, f, M);

  Eigen::Vector2d initsupp = {0.0, 1.0};
  double t = 2.0;

  Eigen::Vector2d speed = {spline.derivative(-1.0), spline.derivative(1.0)};
  Eigen::Vector2d newsupp_ref = initsupp + t * speed;
  Eigen::Vector2d newsupp = findSupport(spline, initsupp, t);

  double error = (newsupp_ref - newsupp).lpNorm<Eigen::Infinity>();
  double tol = 1.e-12;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(CLEmpiricFlux, GodunovFlux_findRoots) {
  double v = -3.0;
  double w = 11.5;
  double z_ref = 1.5;
  auto g = [z_ref](double x) { return std::atan(x - z_ref); };
  double tol = 1.e-4;
  double z = findRoots(v, w, g);
  EXPECT_NEAR(z, z_ref, tol);
}

TEST(CLEmpiricFlux, GodunovFlux_class) {
  Eigen::VectorXd v(10);
  v << -0.14436283, -0.79273996, -0.64006391, -0.61395640, -0.37740092,
      -0.66472778, -0.12799483, -0.08219667, -0.33811969, -0.12848079;

  Eigen::VectorXd w(10);
  w << -0.12204452, -0.21658625, -0.53388862, -0.48357172, -0.03108517,
      -0.33695850, -0.42175670, -0.74430262, -0.20717308, -0.25627699;

  // my solution
  double a = -1.0;
  double b = 0.0;
  unsigned int n = 101;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n + 1, a, b);
  auto f_lambda = [](double x) { return std::sin(PI * x); };
  auto M_lambda = [](double x) { return (-1) * Square(PI) * std::sin(PI * x); };
  Eigen::VectorXd f = x.unaryExpr(f_lambda);
  Eigen::VectorXd M = x.unaryExpr(M_lambda);
  UniformCubicSpline spline(a, b, f, M);
  GodunovFlux godunovFlux(spline);
  Eigen::VectorXd F = v.binaryExpr(w, godunovFlux);

  // reference (computed with code from homework problem
  // "FiniteVolumeSineConsLaw")
  Eigen::VectorXd F_ref(10);
  F_ref << -0.438140682182754, -1.0, -0.994338033898293, -1.0,
      -0.926739695511702, -1.0, -0.39135872034861, -0.25536814568086,
      -0.873445607510785, -0.392763180677777;

  // compare
  double error = (F_ref - F).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-6;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(CLEmpiricFlux, solveCauchyProblem_semiDiscreteRhs) {
  // my solution
  Eigen::VectorXd mu0(10);
  mu0 << 0.0, 0.0, 3.0392, -3.1651, 4.8219, 6.1250, -1.4490, 3.9775, 0.0, 0.0;
  double h = 0.5;
  auto numFlux = [](double v, double w) { return 0.5 * (v + w); };
  Eigen::VectorXd rhs = semiDiscreteRhs(mu0, h, numFlux);

  // reference
  Eigen::VectorXd rhs_ref(10);
  rhs_ref << 0.0, -3.0392, 3.1651, -1.7827, -9.2901, 6.2709, 2.1475, -1.449,
      3.9775, 0.0;

  // compare
  double error = (rhs_ref - rhs).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(CLEmpiricFlux, solveCauchyProblem_RalstonODESolver) {
  // my solution
  Eigen::MatrixXd A(3, 3);
  A << 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
  auto rhs = [A](Eigen::Vector3d y) { return A * y; };
  Eigen::Vector3d y0 = {0.333, -0.75, 0.9};
  double T = PI;
  double n = 200;
  double tau = T / n;
  Eigen::Vector3d yT = RalstonODESolver(rhs, y0, tau, n);

  // reference
  Eigen::MatrixXd expTA(3, 3);
  expTA << std::cos(T), -std::sin(T), 0.0, std::sin(T), std::cos(T), 0.0, 0.0,
      0.0, std::exp(-T);
  Eigen::Vector3d yT_ref = expTA * y0;

  // compare
  double error = (yT_ref - yT).squaredNorm();
  double tol = 1.0e-6;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(CLEmpiricFlux, solveCauchyProblem_computeInitVec) {
  // flux: Burgers' equation
  double a = -1.1;
  double b = 1.1;
  auto f_lambda = [](double u) { return 0.5 * u * u; };
  auto M_lambda = [](double u) { return 1.0; };
  int n = 10;
  Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(n, a, b);
  UniformCubicSpline f(a, b, u.unaryExpr(f_lambda), u.unaryExpr(M_lambda));

  // inital data: jump at zero
  auto u0 = [](double x) { return x < 0.0 ? 1.0 : 0.0; };
  double h = 0.04;
  double T = 1.0;

  // my solution
  Eigen::VectorXd mu0 = computeInitVec(f, u0, h, T);

  // reference
  Eigen::VectorXd mu0_ref(101);
  mu0_ref << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0;

  // compare
  bool same_size = mu0_ref.size() == mu0.size();
  ASSERT_TRUE(same_size);
  if (same_size) {
    double error = (mu0_ref - mu0).lpNorm<Eigen::Infinity>();
    double tol = 1.0e-2;
    EXPECT_NEAR(error, 0.0, tol);
  }
}

TEST(CLEmpiricFlux, solveCauchyProblem_solveCauchyProblem) {
  // flux: Burgers equation
  double a = -1.1;
  double b = 1.1;
  auto f_lambda = [](double u) { return 0.5 * u * u; };
  auto M_lambda = [](double u) { return 1.0; };
  int n = 10;
  Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(n, a, b);
  UniformCubicSpline f(a, b, u.unaryExpr(f_lambda), u.unaryExpr(M_lambda));

  // spacial resolution and final time
  int N = 101;      // number of grid points
  double h = 0.02;  // meshwidth
  double T = 1.0;   // final time

  // inital data
  Eigen::VectorXd mu0(N);
  mu0 << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0;

  // my solution
  Eigen::VectorXd muT = solveCauchyProblem(f, mu0, h, T);

  // reference
  Eigen::VectorXd muT_ref(N);
  muT_ref << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0.999999999999999, 0.999999999999997, 0.999999999999947,
      0.999999999999799, 0.999999999998384, 0.999999999993965,
      0.999999999965757, 0.999999999875356, 0.999999999429691,
      0.999999997971572, 0.999999991902414, 0.999999971744029,
      0.999999896175513, 0.999999642751329, 0.999998746423236,
      0.999995724779585, 0.999985332942479, 0.999950218282813,
      0.999830821891046, 0.999427030793958, 0.998059328803338,
      0.993433061386479, 0.977785682103318, 0.925025173572181,
      0.752433697814842, 0.331614362128096, 0.0224457718791194,
      1.55473915998533e-05, 1.57621528676377e-13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // compare
  double error = (muT_ref - muT).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-4;
  EXPECT_NEAR(error, 0.0, tol);
}

}  // namespace CLEmpiricFlux::test
