/**
 * @file finitevolumesineconslaw_test_mastersolution.cc
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "../finitevolumesineconslaw.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace FiniteVolumeSineConsLaw::test {

TEST(FiniteVolumeSineConsLaw, sineGodFlux) {
  Eigen::VectorXd v(5);
  v << 0.84018771715471, 0.394382926819093, 0.783099223758606,
      0.798440033476073, 0.911647357936784;

  Eigen::VectorXd w(5);
  w << 0.197551369293384, 0.335222755714889, 0.768229594811904,
      0.277774710803188, 0.553969955795431;

  // to test
  Eigen::VectorXd F_GD = v.binaryExpr(w, std::ref(sineGodFlux));

  // reference
  Eigen::VectorXd F_GD_ref(5);
  F_GD_ref << 1.0, 0.945455637559439, 0.665473651823998, 1.0, 0.985660526382152;

  double error = (F_GD - F_GD_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(FiniteVolumeSineConsLaw, sineClawRhs) {
  Eigen::VectorXd v(5);
  v << 0.84018771715471, 0.394382926819093, 0.783099223758606,
      0.798440033476073, 0.911647357936784;

  // to test
  Eigen::VectorXd rhs = FiniteVolumeSineConsLaw::sineClawRhs(v);

  // reference
  Eigen::VectorXd rhs_ref(5);
  rhs_ref << -0.416666666666667, 0.154211752621442, 0.0158953342010465,
      0.132385597521732, -0.30249268434422;

  double error = (rhs - rhs_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(FiniteVolumeSineConsLaw, explTrpzTimestepping) {
  constexpr double LOG2 = 0.693147180559945309417;
  auto g = [](const Eigen::Vector2d& y) { return Eigen::Vector2d(y(1), y(0)); };
  Eigen::Vector2d y0 = {1.0, -1.0};
  unsigned int M = 100;

  // to test
  Eigen::Vector2d sol = explTrpzTimestepping(g, y0, LOG2, M);

  // reference
  auto exact = [](double t) {
    return Eigen::Vector2d(std::exp(-t), -std::exp(-t));
  };
  Eigen::Vector2d sol_ref = exact(LOG2);

  double error = (sol - sol_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-5;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(FiniteVolumeSineConsLaw, solveSineConsLaw) {
  unsigned int N = 10;
  unsigned int M = 10;

  // to test
  Eigen::VectorXd ufinal = solveSineConsLaw(&sineClawRhs, N, M);

  // reference
  Eigen::VectorXd ufinal_ref(10);
  ufinal_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.21656013855545, 0.284395569396954,
      0.22050607120233, 0.14140142198762, 0.0782303913903353;

  double error = (ufinal - ufinal_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(FiniteVolumeSineConsLaw, findTimesteps) {
  unsigned int M = findTimesteps();
  ASSERT_TRUE(155 <= M && M <= 157);
}

TEST(FiniteVolumeSineConsLaw, sineClawReactionRhs) {
  unsigned int N = 10;
  unsigned int M = 10;

  // to test
  double c = 1.0;
  auto bind_c = [c](const Eigen::VectorXd& mu) {
    return sineClawReactionRhs(mu, c);
  };
  Eigen::VectorXd ufinal = solveSineConsLaw(bind_c, N, M);

  // reference
  Eigen::VectorXd ufinal_ref(10);
  ufinal_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0530294982764656, 0.0941150072472929,
      0.0923305241111147, 0.0651004573264405, 0.0366736787456587;

  double error = (ufinal - ufinal_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

}  // namespace FiniteVolumeSineConsLaw::test
