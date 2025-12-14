/**
 * @file pml1d.h
 * @brief NPDE homework 1D Wave Equation with Perfectly Matched Layerscode
 * @author R. Hiptmair
 * @date January 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef PML_H_
#define PML_H_

#define TO_BE_SUPPLEMENTED 1.0

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>

namespace PML1D {
// Default width of PML layer
extern const double L_default;

/** @brief Timestepping for 1D wave equation with PML layers of width L=0.5
 *
 * @param zeta initial state comprising both u and v
 * @param gamma coefficient function gamma sampled on grid
 * @param sigma coefficient function gamma sampled on grid
 * @param M number of timesteps
 * @param T final time
 */
/* SAM_LISTING_BEGIN_1 */
template <typename RECORDER = std::function<void(const Eigen::VectorXd&)>>
Eigen::VectorXd solve1DWavePML(
    const Eigen::VectorXd& zeta_0, const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& sigma, unsigned int M, double T,
    RECORDER&& rec = [](Eigen::VectorXd& /*zeta*/) -> void {}) {
  // Grid resolution parameter N = number of grid nodes - 1
  const unsigned int N = gamma.size() - 1;
  assert(N == sigma.size() - 1);
  assert(2 * N + 1 == zeta_0.size());
  // Vector of basis expansion coefficients, contains
  // $\cob{\vec{\zetabf}^{(k)}}$
  Eigen::VectorXd zeta{zeta_0};
  // Obtain initial velocities at grid nodes by linear interpolation
  Eigen::VectorXd v0 =
      (Eigen::VectorXd(N + 1) << 0.0,
       0.5 * (zeta.segment(N + 1, N - 1) + zeta.segment(N + 2, N - 1)), 0.0)
          .finished();
  // Discretization parameters
  const double L = L_default;          // default width of PML layer
  const double h = (2.0 + 2 * L) / N;  // meshwidth
  const double tau = T / M;            // size of timestep
  // I. Initialize sparse matrices $\cob{\VA,\VR\in\bbR^{2N+1,2N+1}}$.
  Eigen::SparseMatrix<double> A(2 * N + 1, 2 * N + 1);
  Eigen::SparseMatrix<double> R(2 * N + 1, 2 * N + 1);
  // We know that the matrix A has at most 3 non-zero entries per row/column
  A.reserve(Eigen::VectorXi::Constant(2 * N + 1, 3));
  R.reserve(Eigen::VectorXi::Constant(2 * N + 1, 3));
  A.insert(0, 0) = TO_BE_SUPPLEMENTED;
  R.insert(0, 0) = TO_BE_SUPPLEMENTED;
  for (unsigned int i = 1; i < N; ++i) {
    A.insert(i, i) = TO_BE_SUPPLEMENTED;
    R.insert(i, i) = TO_BE_SUPPLEMENTED;
  }
  A.insert(N, N) = TO_BE_SUPPLEMENTED;
  R.insert(N, N) = TO_BE_SUPPLEMENTED;
  for (unsigned int i = 0; i < N; ++i) {
    A.insert(i + N + 1, i + N + 1) = TO_BE_SUPPLEMENTED;
    R.insert(i + N + 1, i + N + 1) = TO_BE_SUPPLEMENTED;
    A.insert(i, i + N + 1) = TO_BE_SUPPLEMENTED;
    A.insert(i + 1, i + N + 1) = TO_BE_SUPPLEMENTED;
    R.insert(i, i + N + 1) = TO_BE_SUPPLEMENTED;
    R.insert(i + 1, i + N + 1) = TO_BE_SUPPLEMENTED;
  }
  for (unsigned int j = 0; j < N; ++j) {
    A.insert(j + N + 1, j) = TO_BE_SUPPLEMENTED;
    A.insert(j + N + 1, j + 1) = TO_BE_SUPPLEMENTED;
    R.insert(j + N + 1, j) = TO_BE_SUPPLEMENTED;
    R.insert(j + N + 1, j + 1) = TO_BE_SUPPLEMENTED;
  }
  // For the sake efficiency Precompute LU factorization of A
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("LU decomposition of A failed");
  }
  // II. Main timestepping loop
  rec(zeta);
  for (unsigned int k = 0; k < M; ++k) {
    // Initialize right-hand side vector $\cob{\vec{\betabf}^{(k)}}$ of size
    // $\cob{2N+1}$
    Eigen::VectorXd beta = -R * zeta;
    beta[0] += 0.5 * h * v0[0];
    for (unsigned int i = 1; i < N; ++i) beta[i] += h * v0[i];
    beta[N] += 0.5 * h * v0[N];
    Eigen::VectorXd zeta_old = zeta;
    zeta = solver.solve(beta);  // $\cob{\vec{\zetabf}^{(k+1)}}$
    rec(zeta);                  // Provide state vector to monitor object
  }
  return zeta;
}
/* SAM_LISTING_END_1 */

/** @brief Energy tracking
 *
 * @param zeta_0 initial vectors of basis expansion coefficients for both u and
 * v
 * @param gamma grid-sampled coefficient
 * @param sigma grid-sampled PML coefficients
 * @param M number of equidistant timesteps
 * @param T final time
 *
 * This function returns the discrete energies for the basis expansion
 * coefficient vectors computed during timestepping.
 */
std::vector<double> trackEnergy(const Eigen::VectorXd& zeta_0,
                                const Eigen::VectorXd& gamma,
                                const Eigen::VectorXd& sigma, unsigned int M,
                                double T);

/** @brief Solving initial value problem and computing errors
 *
 * @param f_v0 functor providing intial v
 * @param f_sol exact solution for u, also used to obtain intial data for u
 * @param N number of spatial grid cells
 * @param M number of timesteps
 * @return approximate L2 and maximum norm errors on [-1,1]
 *
 * @note Uses constant wave speed, quadratic PML profile,
 */
template <typename FUNCTOR_DATA, typename FUNCTOR_SOLUTION>
std::pair<double, double> computeErrorWave1D(FUNCTOR_DATA&& f_v0,
                                             FUNCTOR_SOLUTION&& f_sol,
                                             unsigned int N, unsigned int M) {
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
    const double x = -1.0 - L + i * h;      // Sampling point
    zeta_0[i] = f_sol(x, 0.0);              // Initial value for u
    zeta_0[i + N + 1] = f_v0(x + 0.5 * h);  // Initial value for v
  }
  zeta_0[N] = f_sol(1.0 + L, 0.0);
  // Object for keeping track of approximate solution
  std::vector<Eigen::VectorXd> zetas{};
  auto rec = [&zetas](const Eigen::VectorXd& zeta) -> void {
    zetas.push_back(zeta);
  };
  // Carry out discrete evolution with final time T
  const double T = 4.0;
  const double tau = T / M;
  (void)solve1DWavePML(zeta_0, gamma, sigma, M, T, rec);
  // Compute the error on [-1,1]
  double L2_err = 0.0;
  double Linf_err = 0.0;
  assert(zetas.size() == M + 1);
  Eigen::VectorXd sol(N + 1);  // Contains sampled exact solution
  for (unsigned int k = 1; k <= M; ++k) {
    const double t = k * tau;
    // Sample solution in nodes of the mesh
    for (unsigned int i = 0; i <= N; ++i) {
      const double x = -1.0 - L + i * h;  // Sampling point
      sol[i] = f_sol(x, t);
    }
    // Select index range of nodea lying in the interval [-1,1]
    const int idx_min = static_cast<int>(std::ceil(L / h));
    const int idx_max = static_cast<int>(std::floor((2 + L) / h));
    auto err = sol.segment(idx_min, idx_max - idx_min + 1) -
               zetas[k].segment(idx_min, idx_max - idx_min + 1);
    L2_err += err.squaredNorm();
    Linf_err = std::max(Linf_err, err.lpNorm<Eigen::Infinity>());
  }
  L2_err = std::sqrt(h * tau * L2_err);
  return {L2_err, Linf_err};
}

// Compute and tabulate errors for test case
void tabulateExp1(void);

// Output solution to file for visualization with MATLAB
void plotExp(unsigned int N, unsigned int M, double T,
             std::string filename = "pml1d.m");

}  // namespace PML1D
#endif
