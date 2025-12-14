/**
 * @file fluxlimited_test.cc
 * @brief NPDE homework NewProblem code
 * @author Louis Hurschler
 * @date 26.03.2024
 * @copyright Developed at ETH Zurich
 */

#include "../fluxlimitedfv.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace FluxLimitedFV::test {

TEST(FluxLimitedFV, fluxlimAdvectionTestConstant) {
  const double beta = 1.;
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(100);
  const double h = 0.01;
  const double tau = 0.01;
  const unsigned int nb_timesteps = 100;

  // should stay a vector of 1 because (d mu)/(d x) == 0
  const Eigen::VectorXd mu_res =
      fluxlimAdvection(beta, mu, h, tau, nb_timesteps);
  ASSERT_TRUE(mu_res.size() == 100);
  const Eigen::VectorXd mu_ref = Eigen::VectorXd::Ones(100);

  const double tol = 1.0e-10;
  EXPECT_NEAR(0.0, (mu_res - mu_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FluxLimitedFV, fluxlimAdvectionTestDiscontiuous) {
  const double T = 1.0;
  const double tau = 0.1;
  const double h = 0.1;
  const unsigned N = 10;

  // Parameter coefficient
  const double beta = 1;

  Eigen::VectorXd mu = Eigen::VectorXd::Zero(20);
  for (int i = 0; i < 4; i++) {
    mu(i) = 1.;
  }

  // note that for beta=1. and tau = h, we have (d mu)/(d t) = - (d mu)/(dx)
  // therefore, the discontinuity will simply shift one position to the right in
  // each step
  Eigen::VectorXd flux_sol = fluxlimAdvection(beta, mu, h, tau, N);
  ASSERT_TRUE(flux_sol.size() == 20);

  Eigen::VectorXd ref_sol = Eigen::VectorXd::Zero(20);
  for (int i = 0; i < 14; i++) {
    ref_sol(i) = 1.;
  }

  const double tol = 1.0e-10;
  EXPECT_NEAR(0.0, (flux_sol - ref_sol).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FluxLimitedFV, fluxlimAdvectionTestPhi) {
  const double T = 1.0;
  const double tau = 0.05;
  const double h = 0.05;
  const unsigned N = 20;
  const double beta = 0.8;

  auto phi = [](double theta) {
    return (abs(theta) + theta) / (1.0 + abs(theta));
  };

  Eigen::VectorXd mu(100);
  for (int i = 0; i < 100; i++) {
    double x = i * 0.1;
    if (x < 1) {
      mu(i) = 0.;
    } else if (x < 2) {
      mu(i) = pow(sin(0.5 * M_PI * (x - 1)), 2);
    } else {
      mu(i) = 1.;
    }
  }

  Eigen::VectorXd flux_sol = fluxlimAdvection(beta, mu, h, tau, N, phi);
  ASSERT_TRUE(flux_sol.size() == 100);

  Eigen::VectorXd ref_sol(100);
  ref_sol << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.21816e-30, 6.03115e-27,
      1.87419e-24, 3.67777e-22, 5.11057e-20, 5.34466e-18, 4.36385e-16,
      2.84759e-14, 1.50758e-12, 6.53506e-11, 2.32964e-09, 6.82774e-08,
      1.63391e-06, 3.12622e-05, 0.000453155, 0.00447079, 0.0265037, 0.0915323,
      0.204773, 0.347437, 0.503117, 0.65945, 0.799861, 0.90549, 0.966611,
      0.991738, 0.998663, 0.999873, 0.999995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  const double tol = 1.0e-6;
  EXPECT_NEAR(0.0, (flux_sol - ref_sol).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FluxLimitedFV, fluxlimBurgersTestConstant) {
  Eigen::VectorXd mu0 = Eigen::VectorXd::Ones(100);
  const double h = 0.01;
  const double tau = 0.01;
  const unsigned int nb_timesteps = 100;
  Eigen::VectorXd result = fluxlimBurgers(mu0, h, tau, nb_timesteps);
  ASSERT_TRUE(result.size() == 100);
  const double tol = 1.0e-10;
  ASSERT_NEAR(0.0, (mu0 - result).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FluxLimitedFV, fluxlimBurgersTestPhi) {
  const double h = 0.01;
  const double tau = 0.01;
  const unsigned nb_timesteps = 100;

  auto phi = [](double theta) {
    return (abs(theta) + theta) / (1.0 + abs(theta));
  };

  Eigen::VectorXd mu(100);
  for (int i = 0; i < 100; i++) {
    double x = i * 0.1;
    if (x < 1) {
      mu(i) = 0.;
    } else if (x < 2) {
      mu(i) = pow(sin(0.5 * M_PI * (x - 1)), 2);
    } else {
      mu(i) = 1.;
    }
  }

  Eigen::VectorXd flux_sol = fluxlimBurgers(mu, h, tau, nb_timesteps, phi);
  ASSERT_TRUE(flux_sol.size() == 100);
  Eigen::VectorXd ref_sol(100);
  ref_sol << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00482221, 0.0118288, 0.0196944,
      0.0279691, 0.0364821, 0.0451536, 0.0539392, 0.0628111, 0.0717507,
      0.080745, 0.0897848, 0.0988629, 0.107974, 0.117114, 0.12628, 0.135468,
      0.144676, 0.153904, 0.163148, 0.172407, 0.181681, 0.190967, 0.200266,
      0.209576, 0.218896, 0.228226, 0.237565, 0.246912, 0.256267, 0.265629,
      0.274997, 0.284372, 0.293752, 0.303138, 0.312528, 0.321922, 0.331321,
      0.340723, 0.350128, 0.359536, 0.368947, 0.37836, 0.387775, 0.397191,
      0.406609, 0.416028, 0.425448, 0.434869, 0.444289, 0.45371, 0.46313,
      0.47255, 0.48197, 0.491388, 0.500805, 0.51022, 0.519634, 0.529046,
      0.538456, 0.547863, 0.557268, 0.56667, 0.576068, 0.585463, 0.594855,
      0.604242, 0.613625, 0.623004, 0.632378, 0.641747, 0.65111, 0.660468,
      0.66982, 0.679165, 0.688503, 0.697834, 0.707157, 0.716472, 0.725778,
      0.735075, 0.744363, 0.75364, 0.762905, 0.77216, 0.781401, 0.790629,
      0.799844, 0.80901, 0.819006;

  const double tol = 1.0e-5;
  EXPECT_NEAR(0.0, (flux_sol - ref_sol).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace FluxLimitedFV::test
