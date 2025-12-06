/**
 * @file pml1d.cc
 * @brief NPDE homework 1D Wave Equation with Perfectly Matched Layers code
 * @author R. Hiptmair
 * @date January 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "pml1d.h"

#include <cassert>
#include <fstream>
#include <iomanip>

namespace PML1D {
// The default value for the width of the PML layer
extern const double L_default = 1.0;

void tabulateExp1(void) {
  auto f_u0 = [](double x) -> double {
    if (std::abs(x) > 0.5) {
      return 0.0;
    }
    const double c = std::cos(x * M_PI);
    return c * c;
  };
  auto f_v0 = [](double /*x*/) -> double { return 0.0; };
  auto f_sol = [f_u0](double x, double t) -> double {
    return 0.5 * (f_u0(x + t) + f_u0(x - t));
  };
  // Compute error for different resolutions
  const unsigned int N_min = 8;
  const unsigned int M_min = 8;
  const int L_max = 8;
  // Matrices for storing errors
  Eigen::MatrixXd L2_errs(L_max, L_max);
  Eigen::MatrixXd Linf_errs(L_max, L_max);

  unsigned int N = N_min;
  unsigned int M = M_min;
  std::cout << std::setw(15) << std::fixed << std::setprecision(9) << "N\\M";
  for (unsigned int LM = 0; LM < L_max; ++LM, M *= 2) {
    std::cout << std::setw(15) << M;
  }
  std::cout << std::endl;
  for (unsigned int LN = 0; LN < L_max; ++LN, N *= 2) {
    M = M_min;
    std::cout << std::setw(15) << N << " ";
    for (unsigned int LM = 0; LM < L_max; ++LM, M *= 2) {
      auto [L2_err, Linf_err] = PML1D::computeErrorWave1D(f_v0, f_sol, N, M);
      std::cout << std::setw(15) << L2_err;
      L2_errs(LN, LM) = L2_err;
      Linf_errs(LN, LM) = Linf_err;
    }
    std::cout << std::endl;
  }
}

void plotExp(unsigned int N, unsigned int M, double T, std::string filename) {
  // Exact solution assuming exact absorption
  auto c_gamma = [](double /*x*/) -> double { return 1.0; };
  auto f_u0 = [](double x) -> double {
    if (std::abs(x) > 0.5) {
      return 0.0;
    }
    const double c = std::cos(x * M_PI);
    return c * c;
  };
  auto f_v0 = [](double /*x*/) -> double { return 0.0; };
  auto f_sol = [f_u0](double x, double t) -> double {
    return 0.5 * (f_u0(x + t) + f_u0(x - t));
  };
  // Default PML layer width
  const double L = L_default;
  // PML: coefficient
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
    const double x = -1.0 - L + i * h;      // Sampling point
    zeta_0[i] = f_sol(x, 0.0);              // Initial value for u
    zeta_0[i + N + 1] = f_v0(x + 0.5 * h);  // Initial value for v
  }
  zeta_0[N] = f_sol(1.0 + L, 0.0);
  // Monitor function storing basis expansion coefficient vectors
  std::vector<Eigen::VectorXd> zetas{};
  auto rec = [&zetas, N](const Eigen::VectorXd& zeta) -> void {
    zetas.push_back(zeta.segment(0, N + 1));
  };
  (void)solve1DWavePML(zeta_0, gamma, sigma, M, T, rec);
  // Output solution to file
  std::ofstream outfile(filename);
  std::cout << ">>> Setting: L = " << L << ", N = " << N << ", M = " << M
            << ", T = " << T << "\n zeta_0 = " << zeta_0.transpose() << "\n"
            << std::endl;
  outfile << "data = [ ";
  for (const Eigen::VectorXd& zeta : zetas) {
    outfile << zeta.transpose() << "; ";
  }
  outfile << "];" << std::endl;
}

/* SAM_LISTING_BEGIN_2 */
std::vector<double> trackEnergy(const Eigen::VectorXd& zeta_0,
                                const Eigen::VectorXd& gamma,
                                const Eigen::VectorXd& sigma, unsigned int M,
                                double T) {
  // Grid resolution parameter N = number of grid nodes - 1
  const unsigned int N = gamma.size() - 1;
  const double L = L_default;
  const double h = (2.0 + 2 * L) / N;
  const double tau = T / M;
  // Storage for energies
  std::vector<double> en{};
  // Temporary storage for solution
  Eigen::VectorXd zeta_old(2 * N + 1);
  // Flag telling recorder lambda function that we are at the first timestep
  bool first_step = true;
  // Recorder lambda function
  auto rec = [&](const Eigen::VectorXd& zeta) -> void {
    if (first_step) {
      first_step = false;
    } else {
      auto delta = zeta - zeta_old;
      // Norm of discrete temporal derivative of u-component
      double mu_d_norm = 0.5 * delta[0] * delta[0];
      for (unsigned int i = 1; i < N; ++i) {
        mu_d_norm += delta[i] * delta[i];
      }
      mu_d_norm += 0.5 * delta[N] * delta[N];
      mu_d_norm *= 0.5 * h / (tau * tau);
      // Norm of discrete temporal derivative of v-component
      double nu_d_norm = 0.0;
      for (unsigned int i = 0; i < N; ++i) {
        nu_d_norm += (delta[i + N + 1] * delta[i + N + 1]) / gamma[i];
      }
      nu_d_norm *= 0.5 * h / (tau * tau);
      en.push_back(nu_d_norm + mu_d_norm);
    }
    zeta_old = zeta;
  };
  (void)solve1DWavePML(zeta_0, gamma, sigma, M, T, rec);
  return en;
}
/* SAM_LISTING_END_2 */

}  // namespace PML1D
