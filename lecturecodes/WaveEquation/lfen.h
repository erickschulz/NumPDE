// Demonstration code for course Numerical Methods for Partial Differential
// Equations Author: R. Hiptmair, SAM, ETH Zurich Date: July 2025
// Related to Experiment 9.3.4.20 (Energy conservation for leapfrog)

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#define assertm(exp, msg) assert((void(msg), exp))

namespace LeapfrogWave {

/** @brief Initialization of Poisson matrix
 *
 * Builds the n^2-by-n^2 sparse matrix for the 5‑point Laplacian on an n×n grid
 *  with grid spacing h = 1/(n+1). The matrix A approximates −Δ, i.e.
 * A_{ii} = 4/h^2 and A_{ij} = −1/h^2 for neighbors.
 */
Eigen::SparseMatrix<double> buildLaplacian2D(int n);

/** @brief Initialization of initial vector
 */
template <typename FUNCTOR>
Eigen::VectorXd init_u0(FUNCTOR&& u0_fn, int n) {
  const int N = n * n;       // Total number of vector components
  double h = 1.0 / (n + 1);  // Mesh width
  Eigen::VectorXd u0(N);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const double val = u0_fn((j + 1) * h, (i + 1) * h);
      if (val >= 0.0) {
        u0[(n * i) + j] = val;
      } else {
        u0[(n * i) + j] = 0.0;
      }
    }
  }
  return u0;
}

/** @brief Compute current approximate potential and kinetic energy for leapfrog
 *         timestepping
 *
 */
std::pair<double, double> geten(const Eigen::SparseMatrix<double>& A,
                                double tau, const Eigen::VectorXd& u0,
                                const Eigen::VectorXd& u1);

/** @brief Leapfrof timestepping for constant-coefficients 2D wave equation on
   equidistant spatial mesh: tracking of energies.

   @param u0 initial state (initial velocity set to zero)
   @param n number of grid points in one direction
   @param m number of leapfrog timesteps
 */
std::vector<std::array<double, 4>> leapfrog(const Eigen::VectorXd& u0,
                                            unsigned int n, unsigned int m);

void tabulate_energies(int n, int m, const char* filename = nullptr);

}  // namespace LeapfrogWave
