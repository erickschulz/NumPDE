/**
 * @file contourplot_test.cc
 * @brief NPDE homework ContourPlot code
 * @author Oliver Rietmann
 * @date 25.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../contourplot.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace ContourPlot::test {

constexpr double Square(double x) { return x * x; }

// Leads to isoline F(x) = 2
double F(Eigen::Vector2d x) { return 0.5 * Square(x(0)) + Square(x(1)); };
const Eigen::Vector2d y0(2.0, 0.0);
const double T = 6.0;

// Compute F(x) - 2 along isoline, which should be zero
Eigen::VectorXd errorAlongIsoline(const Eigen::MatrixXd& isolinePoints) {
  int M = isolinePoints.cols();
  Eigen::VectorXd errors(M);
  for (int m = 0; m < M; ++m) {
    Eigen::Vector2d x =
        Eigen::Vector2d(isolinePoints(0, m), isolinePoints(1, m));
    errors(m) = F(x) - 2.0;
  }
  return errors;
}

TEST(ContourPlot, computeIsolinePoints) {
  auto gradF = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(x(0), 2.0 * x(1));
  };
  Eigen::MatrixXd isolinePoints = computeIsolinePoints(gradF, y0, T);
  Eigen::VectorXd errors = errorAlongIsoline(isolinePoints);
  double tol = 1.0e-5;
  ASSERT_NEAR(0.0, errors.lpNorm<Eigen::Infinity>(), tol);
}

TEST(ContourPlot, computeIsolinePointsDQ) {
  Eigen::MatrixXd isolinePointsDQ = computeIsolinePointsDQ(F, y0, T);
  Eigen::VectorXd errors = errorAlongIsoline(isolinePointsDQ);
  double tol = 1.0e-5;
  ASSERT_NEAR(0.0, errors.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace ContourPlot::test
