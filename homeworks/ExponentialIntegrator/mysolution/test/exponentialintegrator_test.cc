/**
 * @file exponentialintegrator_test.cc
 * @brief NPDE homework ExponentialIntegrator
 * @author Tobias Rohner
 * @date 13.04.2021
 * @copyright Developed ate ETH Zurich
 */

#include "../exponentialintegrator.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace ExponentialIntegrator::test {

TEST(exponentialEulerStep, scalar) {
  auto f = [](const Eigen::VectorXd& y) { return y; };
  auto df = [](const Eigen::VectorXd& y) {
    return Eigen::MatrixXd::Ones(1, 1);
  };

  Eigen::VectorXd y0(1);
  y0[0] = 1;
  const double y1 = std::exp(1);

  const Eigen::VectorXd y =
      ExponentialIntegrator::exponentialEulerStep(y0, f, df, 1);

  EXPECT_NEAR(y[0], y1, 1e-8);
}

TEST(exponentialEulerStep, vector) {
  Eigen::MatrixXd A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 0;
  A(1, 0) = 0;
  A(1, 1) = 2;
  auto f = [&](const Eigen::VectorXd& y) { return A * y; };
  auto df = [&](const Eigen::VectorXd& y) { return A; };

  Eigen::VectorXd y0(2);
  y0[0] = 2;
  y0[1] = 1;
  Eigen::VectorXd y1(2);
  y1[0] = std::exp(1) * y0[0];
  y1[1] = std::exp(2) * y0[1];

  const Eigen::VectorXd y =
      ExponentialIntegrator::exponentialEulerStep(y0, f, df, 1);

  EXPECT_NEAR(y[0], y1[0], 1e-8);
  EXPECT_NEAR(y[1], y1[1], 1e-8);
}

}  // namespace ExponentialIntegrator::test
