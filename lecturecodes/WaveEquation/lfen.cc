// Demonstration code for course Numerical Methods for Partial Differential
// Equations Author: R. Hiptmair, SAM, ETH Zurich Date: July 2025
// Related to Experiment 9.3.4.20 (Energy conservation for leapfrog)

#include "lfen.h"

#include <fstream>
#include <iomanip>

namespace LeapfrogWave {

Eigen::SparseMatrix<double> buildLaplacian2D(int n) {
  const int N = n * n;
  Eigen::SparseMatrix<double> A(N, N);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(5 * N);

  const double h = 1.0 / (n + 1);
  const double inv_h2 = 1.0 / (h * h);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int row = (i * n) + j;
      // Diagonal entry
      triplets.emplace_back(row, row, 4.0 * inv_h2);

      // neighbor offsets: up, down, left, right
      if (i > 0) {
        int up = ((i - 1) * n) + j;
        triplets.emplace_back(row, up, -1.0 * inv_h2);
      }
      if (i + 1 < n) {
        int down = ((i + 1) * n) + j;
        triplets.emplace_back(row, down, -1.0 * inv_h2);
      }
      if (j > 0) {
        int left = (i * n) + (j - 1);
        triplets.emplace_back(row, left, -1.0 * inv_h2);
      }
      if (j + 1 < n) {
        int right = (i * n) + (j + 1);
        triplets.emplace_back(row, right, -1.0 * inv_h2);
      }
    }
  }
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

/** @brief Compute current approximate potential and kinetic energy for leapfrog
 *         timestepping
 *
 */
std::pair<double, double> geten(const Eigen::SparseMatrix<double>& A,
                                double tau, const Eigen::VectorXd& u0,
                                const Eigen::VectorXd& u1) {
  auto meanv = 0.5 * (u0 + u1);  // Approximation for velocity
  auto dtemp = (u1 - u0) / tau;  // Approximate temporal derivative
  return {meanv.dot(A * meanv), dtemp.squaredNorm()};
}

/** @brief Leapfrof timestepping for constant-coefficients 2D wave equation on
   equidistant spatial mesh: tracking of energies.

   @param u0 initial state (initial velocity set to zero)
   @param n number of grid points in one direction
   @param m number of leapfrog timesteps
 */
std::vector<std::array<double, 4>> leapfrog(const Eigen::VectorXd& u0,
                                            unsigned int n, unsigned int m) {
  std::vector<std::array<double, 4>> energies;
  const int N = n * n;  // Size of matrices/state vectors
  assertm(u0.size() == N, "Wrong size of initial vector!");
  // Initial velocity
  Eigen::VectorXd v0 = Eigen::VectorXd::Zero(N);

  const double h = 1.0 / (n + 1);  // Mesh width
  // 1. Assemble tridiagonal stiffness matrix as a sparse matrix
  const Eigen::SparseMatrix<double> A = buildLaplacian2D(n);
  // Main timestepping loop: Stormer/Verlet/leapfrog
  double ken;                  // Kinetic energy
  double pen;                  // Potential energy
  const double tau = 1.0 / m;  // Timestep size
  // Special initial step initializes state vector
  Eigen::VectorXd u{u0 + tau * v0 - 0.5 * tau * tau * A * u0};
  std::tie(pen, ken) = geten(A, tau, u0, u);
  energies.push_back(std::array<double, 4>({0.5 * tau, pen, ken, pen + ken}));
  // Auxiliary vectors
  Eigen::VectorXd u_old{u0};
  Eigen::VectorXd u_new(N);
  // Main timestepping loop
  for (int k = 1; k < m; ++k) {
    u_new = -(tau * tau) * (A * u) + 2.0 * u - u_old;
    std::tie(pen, ken) = geten(A, tau, u, u_new);
    energies.push_back(
        std::array<double, 4>({(k + 0.5) * tau, pen, ken, pen + ken}));
    u_old = u;
    u = u_new;
  }
  return energies;
}

void tabulate_energies(int n, int m, const char* filename) {
  auto u0_f = [](double x, double y) -> double {
    const double r =
        std::sqrt(((x - 0.5) * (x - 0.5)) + ((y - 0.5) * (y - 0.5)));
    return (0.2 - r);
  };
  const Eigen::VectorXd u0 = init_u0(u0_f, n);

  auto energies = leapfrog(u0, n, m);
  std::cout << "time t" << std::setw(15) << "pot. en" << std::setw(15)
            << "kin. en" << std::setw(15) << " tot. en." << std::endl;
  std::cout
      << "-----------------------------------------------------------------\n";
  for (const auto& data : energies) {
    std::printf("%1.3f %15f %15f %15f \n", data[0], data[1], data[2], data[3]);
  }
  if (filename) {
    printf("Writing data to %s.csv\n", filename);
    std::ofstream out(std::string(filename) + ".csv");
    if (!out.is_open()) {
      std::cerr << "Error opening file: " << filename << "\n";
      return;
    }
    // Set precision if needed
    out << std::fixed << std::setprecision(16);
    for (const auto& row : energies) {
      out << row[0] << ',' << row[1] << ',' << row[2] << ',' << row[3] << '\n';
    }
    out.close();
  }
}

}  // namespace LeapfrogWave
