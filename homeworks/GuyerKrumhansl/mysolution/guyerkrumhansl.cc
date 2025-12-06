/**
 * @file guyerkrumhansl.cc
 * @brief NPDE homework GuyerKrumhansl code
 * @author R. Hiptmair
 * @date July 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "guyerkrumhansl.h"

#include <iostream>
#include <utility>
#include <vector>

namespace GuyerKrumhansl {

/* Demo code: Triplet initialization of a sparse matrix plus direct solver */
/* SAM_LISTING_BEGIN_1 */
void demo_EigenTripletInit(unsigned int n) {
  const double diag = 2.0;
  const double offd = -1.0;
  // Initialize tridiagonal sparse matrix via COO (Triplet) format
  std::vector<Eigen::Triplet<double>> X_COO;
  X_COO.emplace_back(0, 0, diag);
  X_COO.emplace_back(0, 1, offd);
  for (int j = 1; j < n - 1; ++j) {
    X_COO.emplace_back(j, j - 1, offd);
    X_COO.emplace_back(j, j, diag);
    X_COO.emplace_back(j, j + 1, offd);
  }
  X_COO.emplace_back(n - 1, n - 1, diag);
  X_COO.emplace_back(n - 1, n - 2, offd);
  // Build sparse matrix
  Eigen::SparseMatrix<double> X(n, n);
  X.setFromTriplets(X_COO.begin(), X_COO.end());
  // Solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(X);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not factorize X");
  }
  std::cout << "x = "
            << solver.solve(Eigen::VectorXd::Constant(n, 1.0)).transpose()
            << std::endl;
}
/* SAM_LISTING_END_1 */

/* Demo code: efficient direct initialization of a sparse matrix */
/* SAM_LISTING_BEGIN_9 */
void demo_EigenSparseInit(unsigned int n, unsigned int m) {
  // Set up zero sparse matrix with row major storage format
  // This format is essential for being able to set the maximal
  // number of non-zero entries \textbf{per row}.
  Eigen::SparseMatrix<double, Eigen::RowMajor> X(n, m);
  // Reserve space for at most nnz\_row non-zero entries per row
  const std::size_t nnz_row = 3;
  X.reserve(Eigen::VectorXi::Constant(n, nnz_row));
  // Initialize nnz\_row  entries per row
  for (int row_idx = 0; row_idx < n; ++row_idx) {
    for (int k = 0; k < nnz_row; ++k) {
      X.insert(row_idx, (row_idx * k) % m) = 1.0;
    }
  }
  Eigen::MatrixXd X_dense = X;
  std::cout << "Matrix X = " << std::endl << X_dense << std::endl;
}
/* SAM_LISTING_END_1 */

#define MISSING_VALUE (-42.0)

// Initialization of matrices for MOL ODE
/* SAM_LISTING_BEGIN_2 */
std::pair<std::vector<Eigen::Triplet<double>>,
          std::vector<Eigen::Triplet<double>>>
generateMOLMatrices(unsigned int M, double rho, double sigma, double mu,
                    double kappa) {
  const unsigned int N = 2 * M - 1;  // Total number of d.o.f.s
  const double h = 1.0 / M;          // Meshwidth
  // Triplet vectors to be returned
  std::vector<Eigen::Triplet<double>> Mh_COO;
  std::vector<Eigen::Triplet<double>> Ah_COO;
// I. Matrix $\wh{\VM}$
  for (int j = 0; j < M; ++j) {
    Mh_COO.emplace_back(j, j, rho * h);
  }
  const double fac = sigma * h / 6.0;
  Mh_COO.emplace_back(M, M, 4.0 * fac);
  Mh_COO.emplace_back(M, M + 1, fac);
  for (int j = 1; j < M - 2; ++j) {
    Mh_COO.emplace_back(j + M, j + M, 4 * fac);
    Mh_COO.emplace_back(j + M, j + M - 1, fac);
    Mh_COO.emplace_back(j + M, j + M + 1, fac);
  }
  Mh_COO.emplace_back(N - 1, N - 1, 4 * fac);
  Mh_COO.emplace_back(N - 1, N - 2, fac);
// II. Matrix $\wh{\VA}$
  Ah_COO.emplace_back(0, M, -1.0);
  Ah_COO.emplace_back(M, 0, kappa);
  for (int j = 1; j < M - 1; ++j) {
// ************************
// Your code here
// ************************
  }
  Ah_COO.emplace_back(M - 1, N - 1, 1.0);
  Ah_COO.emplace_back(N - 1, M - 1, -kappa);
  const double diag =
      // ************************
      // Replace this line
      MISSING_VALUE;
  // ************************
  const double offd =
      // ************************
      // Replace this line
      MISSING_VALUE;
  // ************************
  Ah_COO.emplace_back(M, M, diag);
  Ah_COO.emplace_back(M, M + 1, offd);
  for (int j = 1; j < M - 2; ++j) {
// ************************
// Your code here
// ************************
  }
  Ah_COO.emplace_back(N - 1, N - 1, diag);
  Ah_COO.emplace_back(N - 1, N - 2, offd);

  return {Mh_COO, Ah_COO};
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::vector<double> track_GKHeatEnergy(double T, unsigned int L, unsigned int M,
                                       const Eigen::VectorXd& mu0,
                                       const Eigen::VectorXd& zeta0, double rho,
                                       double sigma, double mu, double kappa) {
  std::vector<double> energies;
  // ************************
  // Extend to a meaningful implementation
  auto rec = [/* Your captured local variables here */](
                 const Eigen::VectorXd& nu) -> void { /* Your code here*/ };
// ************************
  (void)timestepping_GKHeat(T, L, M, mu0, zeta0, rho, sigma, mu, kappa, rec);
  return energies;
}
/* SAM_LISTING_END_3 */

}  // namespace GuyerKrumhansl
