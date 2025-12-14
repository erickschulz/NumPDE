/**
 * @file pml1d_test.cc
 * @brief NPDE homework PML1D code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../pml1d.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace PML1D::test {

static constexpr double L_default = 1;

// Analytic solution
double f_sol(double x, double t) {
  auto f_u0 = [](double x) -> double {
    if (std::abs(x) > 0.5) {
      return 0.0;
    }
    const double c = std::cos(x * M_PI);
    return c * c;
  };
  return 0.5 * (f_u0(x + t) + f_u0(x - t));
}

std::vector<Eigen::VectorXd> solve(unsigned N, unsigned M, double T) {
  // Analytic solution
  // PML coefficient function
  const double s0 = 10.0;
  auto c_sigma = [s0](double x) -> double {
    if (std::abs(x) <= 1.0) {
      return 0.0;
    }
    if (x < -1.0) {
      return s0 * (x + 1.0) * (x + 1.0);
    }
    return s0 * (x - 1.0) * (x - 1.0);
  };
  // Wavespeed coefficient function
  auto c_gamma = [](double /*x*/) -> double { return 1.0; };
  // Default PML layer width
  const double L = L_default;
  // Sampling coefficients on equidistant mesh
  const double h = (2.0 + 2 * L) / N;
  Eigen::VectorXd sigma(N + 1);
  Eigen::VectorXd gamma(N + 1);
  for (unsigned int i = 0; i <= N; ++i) {
    const double x = -1.0 - L + i * h;  // Sampling point
    sigma[i] = c_sigma(x);
    gamma[i] = c_gamma(x);
  }
  // Sampling initial data on equidistant mesh
  Eigen::VectorXd zeta_0(2 * N + 1);
  for (unsigned int i = 0; i < N; ++i) {
    const double x = -1.0 - L + i * h;  // Sampling point
    zeta_0[i] = f_sol(x, 0.0);          // Initial value for u
    zeta_0[i + N + 1] = 0;              // Initial value for v
  }
  zeta_0[N] = f_sol(1.0 + L, 0.0);
  // Object for keeping track of approximate solution
  std::vector<Eigen::VectorXd> zetas;
  auto rec = [&](const Eigen::VectorXd& zeta) -> void {
    zetas.push_back(zeta);
  };
  // Carry out discrete evolution with final time T
  solve1DWavePML(zeta_0, gamma, sigma, M, T, rec);
  return zetas;
}

TEST(PML1D, wave) {
  const double L = L_default;
  const unsigned N = 512;
  const double h = (2 * L + 2) / N;
  const unsigned M = N;
  const double T = 0.5;
  const double tau = T / M;
  const auto zetas = solve(N, M, T);
  Eigen::VectorXd ref_sol(N + 1);  // Contains sampled exact solution
  double Linf_err = 0;
  for (unsigned int k = 1; k <= M; ++k) {
    const double t = k * tau;
    // Sample solution in nodes of the mesh
    for (unsigned int i = 0; i <= N; ++i) {
      const double x = -1.0 - L + i * h;  // Sampling point
      ref_sol[i] = f_sol(x, t);
    }
    // Select index range of nodea lying in the interval [-1,1]
    const int idx_min = static_cast<int>(std::ceil(L / h));
    const int idx_max = static_cast<int>(std::floor((2 + L) / h));
    auto err = ref_sol.segment(idx_min, idx_max - idx_min + 1) -
               zetas[k].segment(idx_min, idx_max - idx_min + 1);
    Linf_err = std::max(Linf_err, err.lpNorm<Eigen::Infinity>());
  }
  ASSERT_NEAR(Linf_err, 0, 1e-3);
}

TEST(PML1D, boundarynodes) {
  const unsigned N = 512;
  const unsigned M = N;
  const double T = 4;
  const auto zetas = solve(N, M, T);
  for (unsigned k = 1; k <= M; ++k) {
    ASSERT_NEAR(zetas[k][0], 0, 1e-1);
    ASSERT_NEAR(zetas[k][N], 0, 1e-1);
  }
}

}  // namespace PML1D::test
