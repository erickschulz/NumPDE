/**
 * @ file LinearFE1D.h
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher, Erick Schulz
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

#include <Eigen/SparseCholesky>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace LinearFE1D {

// Calculate the matrix entries corresponding to the Galerkin matrix
// of the Laplace operator with non-constant coefficient function alpha
// using composite midpoint integration rule. The entries are stored in
// the Eigen triplet data structure (COO sparse matrix format).
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeA(const Eigen::VectorXd& mesh,
                                             FUNCTOR1&& alpha) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializing the vector of triplets whose size corresponds to the
  // number of entries in a (N+1) x (N+1) tridiagonal (band) matrix
  std::vector<Eigen::Triplet<double>> triplets;
  // Maximal size of vector; avoids intermediate allocations
  triplets.reserve(3 * N + 1);
  //====================
  // Your code goes here
  //====================
  return triplets;
}  // computeA
/* SAM_LISTING_END_1 */

// Calculate the matrix entries corresponding to the mass matrix
// associated to the L2 pairing weighted by the variable coefficient
// function gamma using the trapezoidal integration rule.
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeM(const Eigen::VectorXd& mesh,
                                             FUNCTOR1&& gamma) {
  // Nodes are indexed as $0=x_0 < x_1 < ... < x_N = 1$
  unsigned N = mesh.size() - 1;

  // The basis tent functions spanning the space of continuous piecewise
  // linear polynomial satisfy the cardinal basis property with respect
  // to the nodes of the mesh. Using  the trapezoidal rule to approximate
  // the integrals therefore lead to a diagonal $(N+1) \times (N+1)$ matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N + 1);

  //====================
  // Your code goes here
  //====================

  return triplets;
}  // computeM
/* SAM_LISTING_END_2 */

// Calculate the entries of the right hand side vector using
// the composite trapezoidal integration rule
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR1>
Eigen::VectorXd computeRHS(const Eigen::VectorXd& mesh, FUNCTOR1&& f) {
  // Nodes are indexed as $0=x_0 < x_1 < ... < x_N = 1$
  unsigned N = mesh.size() - 1;
  // Initializing right hand side vector
  Eigen::VectorXd rhs_vec = Eigen::VectorXd::Zero(N + 1);

  //====================
  // Your code goes here
  //====================

  return rhs_vec;
}  // computeRHS
/* SAM_LISTING_END_3 */

// SOLVE THE LINEAR SYSTEM OF PROBLEM (A)
/* SAM_LISTING_BEGIN_A */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveA(const Eigen::VectorXd& mesh, FUNCTOR1&& gamma,
                       FUNCTOR2&& f) {
  // Nodes are indexed as $0=x_0 < x_1 < ... < x_N = 1$
  // Note: What is called N here, is M in the project description!
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u = Eigen::VectorXd::Zero(N + 1);  // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);       // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N + 1, N + 1);       // mass galerkin mat
  Eigen::SparseMatrix<double> L(N + 1, N + 1);       // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce (zero) dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE L*u = rhs\_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  // The solution vector u was initialized with zeros, and therefore already
  // contains the zero Dirichlet boundary data in the first and last entry
  return u;
}  // solveA
/* SAM_LISTING_END_A */

// SOLVE THE LINEAR SYSTEM OF PROBLEM (B)
/* SAM_LISTING_BEGIN_B */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveB(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& f, double u0, double u1) {
  // Nodes are indexed as $0=x_0 < x_1 < ... < x_N = 1$
  unsigned N = mesh.size() - 1;
  // Initializations
  Eigen::VectorXd u(N + 1);                     // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);  // laplacian galerkin mat
  // Some tool variables
  double dx_left, dx_right;  // cell widths

  // I. Build the Galerkin matrix A
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs\_vec using an built-in Eigen solver
  //====================
  // Your code goes here
  //====================

  // The solution vector still needs to be supplemented with the known
  // boundary values
  u(0) = u0;  // left boundary node
  u(N) = u1;  // right boundary node
  return u;
}  // solveB
/* SAM_LISTING_END_B */

// Build an sol!ve the LSE corresponding to (C)
/* SAM_LISTING_BEGIN_C */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveC(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& gamma) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u(N + 1);                     // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);  // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N + 1, N + 1);  // mass galerkin mat
  Eigen::SparseMatrix<double> L(N + 1, N + 1);  // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs\_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  return u;
}  // solveC
/* SAM_LISTING_END_C */

}  // namespace LinearFE1D
