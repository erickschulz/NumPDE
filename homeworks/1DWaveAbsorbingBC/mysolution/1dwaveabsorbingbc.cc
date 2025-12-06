/**
 * @file 1dwaveabsorbingbc.cc
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "1dwaveabsorbingbc.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <utility>
#include <vector>

namespace WaveAbsorbingBC1D {

constexpr double PI = 3.14159265358979323846;

// Boundary data on the right-hand side (x = 1)
double g(double t) { return 0 <= t && t <= PI ? std::sin(t) : 0.0; }

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix A
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param c speed of propagation
 * @param h (spacial) meshwidth
 * @return Full Galerkin matrix A of shape (N+1)x(N+1)
 */
/* SAM_LISTING_BEGIN_7 */
Eigen::SparseMatrix<double> getA_full(unsigned int N, double c, double h) {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * (N + 1) - 2);  // that many triplets needed
  const double scale = c * c / h;
  // store first row separately
  triplets.push_back(Eigen::Triplet<double>(0, 0, scale));
  triplets.push_back(Eigen::Triplet<double>(0, 1, -scale));
  // loop over interior rows
  for (unsigned i = 1; i < N; ++i) {
    triplets.push_back(Eigen::Triplet<double>(i, i - 1, -scale));
    triplets.push_back(Eigen::Triplet<double>(i, i, 2.0 * scale));
    triplets.push_back(Eigen::Triplet<double>(i, i + 1, -scale));
  }
  // store last row separately
  triplets.push_back(Eigen::Triplet<double>(N, N - 1, -scale));
  triplets.push_back(Eigen::Triplet<double>(N, N, scale));
  // Creat (N+1)x(N+1) sparse matrix in CRS format
  Eigen::SparseMatrix<double> A(N + 1, N + 1);
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}
/* SAM_LISTING_END_7 */

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix B
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param c speed of propagation
 * @return Full Galerkin matrix B of size (N+1)x(N+1)
 */
/* SAM_LISTING_BEGIN_8 */
Eigen::SparseMatrix<double> getB_full(unsigned int N, double c) {
  Eigen::SparseMatrix<double> B(N + 1, N + 1);
  // Just a single non-zero entry; we can afford to sete it directly
  B.coeffRef(0, 0) = c;
  return B;
}
/* SAM_LISTING_END_8 */

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix M
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param h (spacial) meshwidth
 * @return Full Galerkin matrix M of size (N+1)x(N+1)
 * Note that M is a diagonal matrix!
 */
/* SAM_LISTING_BEGIN_9 */
Eigen::SparseMatrix<double> getM_full(unsigned int N, double h) {
  Eigen::SparseMatrix<double> M(N + 1, N + 1);
  M.setIdentity();  // Supposed to be efficient
  M *= h;
  // Modify two entries; efficiency does not matter much
  M.coeffRef(0, 0) = h / 2;
  M.coeffRef(N, N) = h / 2;
  return M;
}
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd waveLeapfrogABC(double c, double T, unsigned int N,
                                unsigned int m) {
  // N is also the number of cells of the mesh
  double h = 1.0 / N;
  // Obtain Galerkin matrices for truncated finite element space
  // Note that the functions get*_full return the matrices for the full finite
  // element space including the tent function located at x=1. Removing the last
  // row and column of that matrix amounts to dropping that basis function.
  // However, the efficiency of this block() operation in the case of sparse
  // matrices is in doubt, in particular, since the result is assigned to
  // another sparse matrix, which foils Eigen's expression template
  // optimization. The use of "auto" would be highly advisable here!
  Eigen::SparseMatrix<double> A = getA_full(N, c, h).block(0, 0, N, N);
  Eigen::SparseMatrix<double> B = getB_full(N, c).block(0, 0, N, N);
  Eigen::SparseMatrix<double> M = getM_full(N, h).block(0, 0, N, N);
  // Matrix for returning solution
  Eigen::MatrixXd R(m + 1, N + 1);
  //====================
  // Your code goes here
  //====================
  return R;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> computeEnergies(
    const Eigen::MatrixXd& full_solution, double c, double tau) {
  int m = full_solution.rows() - 1;
  int N = full_solution.cols() - 1;
  double h = 1.0 / N;

  Eigen::SparseMatrix<double> A = getA_full(N, c, h);
  Eigen::SparseMatrix<double> M = getM_full(N, h);

  Eigen::VectorXd E_pot(m + 1);
  Eigen::VectorXd E_kin(m);

  //====================
  // Your code goes here
  //====================

  return std::make_pair(E_pot, E_kin);
}
/* SAM_LISTING_END_2 */

}  // namespace WaveAbsorbingBC1D
