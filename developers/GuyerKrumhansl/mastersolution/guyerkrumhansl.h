/**
 * @file guyerkrumhansl.h
 * @brief NPDE homework GuyerKrumhansl code
 * @author R. Hiptmair
 * @date July 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cassert>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace GuyerKrumhansl {

/** @brief Demo code: Triplet initialization of a sparse matrix plus direct
 * solver */
void demo_EigenTripletInit(unsigned int n);

/** @brief Demo code: efficient direct initialization of a sparse matrix
    @param m,n: matrix dimensions
 */
void demo_EigenSparseInit(unsigned int n = 20, unsigned int m = 10);

/** @brief Building of the matrices for the MOL ODE */
std::pair<std::vector<Eigen::Triplet<double>>,
          std::vector<Eigen::Triplet<double>>>
generateMOLMatrices(unsigned int M, double rho, double sigma, double mu,
                    double kappa);

/** @brief Implicit Euler timestepping for Guyer-Krumhansl heat equation */
/* SAM_LISTING_BEGIN_1 */
template <typename RECORDER = std::function<void(const Eigen::VectorXd&)>>
std::pair<Eigen::VectorXd, Eigen::VectorXd> timestepping_GKHeat(
    double T, unsigned int L, unsigned int M, const Eigen::VectorXd& mu0,
    const Eigen::VectorXd& zeta0, double rho, double sigma, double mu,
    double kappa, RECORDER&& rec = [](const Eigen::VectorXd& /*nu*/) {}) {
  const double tau = T / L;          // size of timestep
  const unsigned int N = 2 * M - 1;  // Total number of FE d.o.f.s

  assert((void("Size mismatch for mu0"), mu0.size() == M));
  assert((void("Size mismatch for zeta0"), zeta0.size() == M - 1));
  // D.o.f. vector for timestepping
  Eigen::VectorXd nu(N);
  nu << mu0, zeta0;
  (void)rec(nu);

  // Obtain triplet vectors for MOL ODE matrices
  auto [Mh_COO, Ah_COO] = generateMOLMatrices(M, rho, sigma, mu, kappa);
  Eigen::SparseMatrix<double> Mh(N, N);
  Eigen::SparseMatrix<double> Xh(N, N);
#if SOLUTION
  // Fill sparse matrix $\wh{\VM}$
  Mh.setFromTriplets(Mh_COO.begin(), Mh_COO.end());
  // Build matrix $\wh{\VX} := \wh{\VM}-\tau \wh{\VA}$
  std::vector<Eigen::Triplet<double>> Xh_COO = Mh_COO;
  for (Eigen::Triplet<double>& t : Ah_COO) {
    Xh_COO.emplace_back(t.row(), t.col(), -tau * t.value());
  }
  // Build sparse matrix in CRS format
  Xh.setFromTriplets(Xh_COO.begin(), Xh_COO.end());
#else
  // ************************
  // Your code here
  // ************************
#endif
  // LU factorization of $\wh{\VX}$
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(Xh);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not factorize Xh");
  }
  // Main timestepping loop
  for (int k = 0; k < L; ++k) {
    // Implement \prbeqref{eq:eerec}
    nu = solver.solve(Mh * nu);
    (void)rec(nu);
  }
  return {nu.head(M), nu.tail(M - 1)};
}
/* SAM_LISTING_END_1 */

/** @brief tracking evolution of discrete energy during timestepping */
std::vector<double> track_GKHeatEnergy(double T, unsigned int L, unsigned int M,
                                       const Eigen::VectorXd& mu0,
                                       const Eigen::VectorXd& zeta0, double rho,
                                       double sigma, double mu, double kappa);

}  // namespace GuyerKrumhansl
