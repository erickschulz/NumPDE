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
#include <fstream>

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
  u_ref << -4.110139126548574e-18, -4.726609052952856e-18,
      8.005812256622144e-18, 4.112943916130766e-18, -2.229801237005196e-18,
      -1.058084400587765e-17, -5.381464657340299e-18, -5.003786039870063e-18,
      2.673653225362404e-18, 2.840941679419017e-18, -1.252352592776204e-18,
      2.490809359404724e-18, -4.201777377022271e-18, 5.706155480941125e-18,
      5.068577032592605e-18, -9.153078596435666e-19, -1.221190297957524e-18,
      7.468187804938338e-19, -5.402764037935826e-18, 3.077497252877032e-18,
      4.418893015442913e-18, -2.219883920212623e-18, 1.570773367582107e-18,
      -4.442373871064228e-18, 8.620699240530914e-19, 9.259701980933254e-18,
      -1.031351217221853e-17, 1.655713505636003e-18, 4.066975173414383e-18,
      -8.248640558190562e-18, 4.522462362553913e-18, 5.024357001909054e-18,
      -4.807505389136438e-18, 9.159392103825692e-18, -1.163367305239895e-17,
      4.261429897291949e-18, -1.717464766284146e-18, 1.344575824044138e-18,
      -6.790372340157351e-19, 0.02195039711368491, 0.995324251295099,
      0.9722962073901446, 0.9508328648027865, 0.9284662194065862,
      0.9055773116931443, 0.8823479253068323, 0.858875149626509,
      0.8352157687675517, 0.8114051618558736, 0.7874661908282345,
      0.7634137333012438, 0.7392571090124367, 0.7150013648374348,
      0.6906478416314827, 0.6661941705428048, 0.6416336391642492,
      0.6169535780357457, 0.5921317368829708, 0.5671273872983082,
      0.5418535842941374, 0.5160284069127383, 0.4839715930872616,
      0.4581464157058627, 0.4328726127016917, 0.4078682631170292,
      0.3830464219642545, 0.3583663608357506, 0.3338058294571951,
      0.3093521583685174, 0.2849986351625652, 0.2607428909875631,
      0.236586266698756, 0.212533809171765, 0.1885948381441263,
      0.1647842312324484, 0.1411248503734911, 0.1176520746931679,
      0.09442268830685582, 0.07153378059341382, 0.04916713519721346,
      0.02770379260985541;

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
