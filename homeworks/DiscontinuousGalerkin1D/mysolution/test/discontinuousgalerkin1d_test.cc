/**
 * @file discontinuousgalerkin1d_test_mastersolution.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include "../discontinuousgalerkin1d.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace DiscontinuousGalerkin1D::test {

TEST(DiscontinuousGalerkin1D, compBmat) {
  int Ml = 1;
  int Mr = 2;
  double h = 0.55;

  // to test
  Eigen::MatrixXd B(compBmat(Ml, Mr, h));

  // reference
  int N = 2 * (Ml + Mr + 1);
  Eigen::VectorXd d(N);
  for (int i = 0; i < N; i += 2) {
    d[i] = 0.55;
    d[i + 1] = 0.013864583333333;
  }
  Eigen::MatrixXd B_ref = d.asDiagonal();

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (B - B_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, G) {
  Eigen::VectorXd mu(6);
  mu << 0.1, 0.2, 0.3, 0.4, 0.3, 0.2;
  auto f = [](double x) { return 0.5 * x * x; };
  int Ml = 1;
  int Mr = 1;
  double h = 0.55;
  double tol = 1.0e-8;

  Eigen::VectorXd Gvec;
  Eigen::VectorXd Gvec_ref(6);

  auto Fv = [](double v, double w) { return v; };
  Gvec = G(mu, f, Fv, Ml, Mr, h);
  Gvec_ref << 0.155, 0.0395977083333333, 0.255, 0.129515833333333, -0.055,
      0.185347708333333;
  ASSERT_NEAR(0.0, (Gvec - Gvec_ref).lpNorm<Eigen::Infinity>(), tol);

  auto Fw = [](double v, double w) { return w; };
  Gvec = G(mu, f, Fw, Ml, Mr, h);
  Gvec_ref << 0.145, 0.0615977083333333, 0.055, 0.0937658333333333, -0.245,
      0.0423477083333333;
  ASSERT_NEAR(0.0, (Gvec - Gvec_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, dgcl) {
  Eigen::VectorXd mu0(6);
  mu0 << 0.1, 0.2, 0.3, 0.4, 0.3, 0.2;
  auto f = [](double x) { return 0.5 * x * x; };
  auto F = [](double v, double w) { return v; };
  int Ml = 1;
  int Mr = 1;
  double h = 0.55;
  double T = 1.0;
  unsigned int m = 2;

  Eigen::VectorXd mu =
      DiscontinuousGalerkin1D::dgcl(mu0, f, F, T, Ml, Mr, h, m);

  Eigen::VectorXd mu_ref(6);
  mu_ref << 0.524398178802972, 5.22582383256893, 0.788213539122826,
      17.3887446322133, -0.540429592450715, 17.1803412500935;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu - mu_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, Feo) {
  Eigen::Vector4d v = {1.0, -1.0, 1.0, 0.3};
  Eigen::Vector4d w = {2.0, 2.0, 1.0, -2.0};

  Eigen::Vector4d F = v.binaryExpr(w, std::ref(Feo));
  Eigen::Vector4d F_ref = {-2.0, -4.25, 0.0, 0.21};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (F - F_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, solveTrafficFlow) {
  int Ml = 40;
  int Mr = 40;
  int N_half = Mr + Ml + 1;

  Eigen::VectorXd x_ref = Eigen::VectorXd::LinSpaced(N_half, -2.0, 2.0);

  Eigen::VectorXd u_ref(N_half);
  u_ref << 3.64755156791648e-18, 4.39830640427375e-18, -2.31796116971711e-18, 2.9538987603557e-18, -1.58521014805857e-19, 3.20751082359005e-18, -6.4724349998136e-18, -2.30483923279167e-18, 1.39383490232115e-18, -3.22882041631168e-18, 1.08384064850155e-18, -1.45610350430611e-19, 2.6336678175536e-18, -1.39750028627349e-18, -1.38788841774047e-18, -3.1418630767741e-18, -1.37658662527077e-18, -1.85677200227667e-18, 9.59174622147387e-18, -7.89294336958034e-18, 1.16799052452298e-19, 5.03244526876313e-19, 7.55130197778723e-19, 7.8800859632695e-18, -6.09424981313119e-18, -1.42256976787597e-18, 9.37099209110105e-18, -1.10054805031767e-17, -3.12086671631005e-22, 1.7173522027843e-18, 2.90373304148436e-18, -3.34680065653073e-18, 7.36365255110097e-18, -1.77977892422898e-18, -9.8298778151396e-21, -2.33375638543984e-18, 6.22334136481836e-18, -8.00181107777641e-18, -4.37230711530437e-18, 0.0177040238183212, 0.971866831980608, 0.950832864802787, 0.928466219406586, 0.905577311693144, 0.882347925306832, 0.858875149626509, 0.835215768767552, 0.811405161855874, 0.787466190828235, 0.763413733301244, 0.739257109012437, 0.715001364837435, 0.690647841631483, 0.666194170542805, 0.641633639164249, 0.616953578035746, 0.592131736882971, 0.567127387298308, 0.541853584294138, 0.516028406912739, 0.483971593087262, 0.458146415705863, 0.432872612701692, 0.407868263117029, 0.383046421964254, 0.35836636083575, 0.333805829457195, 0.309352158368517, 0.284998635162565, 0.260742890987563, 0.236586266698756, 0.212533809171765, 0.188594838144126, 0.164784231232449, 0.141124850373491, 0.117652074693168, 0.0944226883068558, 0.0715337805934138, 0.0491671351972134, 0.0277037926098554, 0.0080342111763938;

  DiscontinuousGalerkin1D::Solution solution = solveTrafficFlow();

  // test for correct length
  ASSERT_EQ(solution.x_.size(), solution.u_.size());
  ASSERT_EQ(solution.x_.size(), N_half);

  // test for correct values
  if (solution.x_.size() == N_half && solution.u_.size() == N_half) {
    double tol = 1.0e-5;
    ASSERT_NEAR(0.0, (solution.x_ - x_ref).lpNorm<Eigen::Infinity>(), tol);
    ASSERT_NEAR(0.0, (solution.u_ - u_ref).lpNorm<Eigen::Infinity>(), tol);
  }
}

}  // namespace DiscontinuousGalerkin1D::test
