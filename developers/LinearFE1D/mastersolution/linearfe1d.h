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
#if SOLUTION
  // Some auxiliary variables
  double diag, off_diag;     // matrix entries
  double dx_left, dx_right;  // cell widths

  /* Computing diagonal entries */
  // First diagonal entry (left boundary node)
  dx_right = mesh(1) - mesh(0);
  diag = alpha((mesh(1) + mesh(0)) / 2.) / dx_right;
  triplets.push_back(Eigen::Triplet<double>(0, 0, diag));
  // Last diagonal entry (right boundary node)
  dx_left = mesh(N) - mesh(N - 1);
  diag = alpha((mesh(N) + mesh(N - 1)) / 2.) / dx_left;
  triplets.push_back(Eigen::Triplet<double>(N, N, diag));
  // All diagonal entries associated to interior nodes
  for (unsigned i = 1; i < N; ++i) {
    dx_left = mesh(i) - mesh(i - 1);
    dx_right = mesh(i + 1) - mesh(i);
    diag = alpha((mesh(i - 1) + mesh(i)) / 2.) / dx_left +
           alpha((mesh(i) + mesh(i + 1)) / 2.) / dx_right;
    triplets.push_back(Eigen::Triplet<double>(i, i, diag));
  }

  /* Computing off-diagonal entries */
  for (unsigned i = 0; i < N; ++i) {
    dx_right = mesh(i + 1) - mesh(i);
    off_diag = -alpha((mesh(i) + mesh(i + 1)) / 2.) / dx_right;
    triplets.push_back(Eigen::Triplet<double>(i + 1, i, off_diag));
    triplets.push_back(Eigen::Triplet<double>(i, i + 1, off_diag));
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
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

#if SOLUTION
  // Some tool variables
  double diag, off_diag;
  double dx;  // cell widths

  /* Computing diagonal entries */
  // First diagonal entry (left boundary node)
  dx = mesh(1) - mesh(0);
  diag = gamma(mesh(0)) * 0.5 * dx;
  triplets.push_back(Eigen::Triplet<double>(0, 0, diag));
  // Last diagonal entry (right boundary node)
  dx = mesh(N) - mesh(N - 1);
  diag = gamma(mesh(N)) * 0.5 * dx;
  triplets.push_back(Eigen::Triplet<double>(N, N, diag));
  // All diagonal entries associated to interior nodes
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    diag = gamma(mesh(i)) * 0.5 * dx;
    triplets.push_back(Eigen::Triplet<double>(i, i, diag));
  }  // computeM
#else
  //====================
  // Your code goes here
  //====================
#endif

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

#if SOLUTION
  // Some tool variables
  double dx;

  /* Calculate the entries */
  // First entry (left boundary node)
  rhs_vec(0) = f(mesh(0)) * (mesh(1) - mesh(0)) * 0.5;
  // Last entry (right boundary node)
  rhs_vec(N) = f(mesh(N)) * (mesh(N) - mesh(N - 1)) * 0.5;
  // All other enties associated to interior nodes
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    rhs_vec(i) = f(mesh(i)) * 0.5 * dx;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

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
#if SOLUTION
  // I.i Compute the entries of the Laplace Galerkin matrix A
  // (const. ceoff. func. alpha = 1.0)
  std::vector<Eigen::Triplet<double>> triplets_A =
      computeA(mesh, [](double) -> double { return 1.0; });
  // I.ii Compute the entries of the mass matrix M
  std::vector<Eigen::Triplet<double>> triplets_M = computeM(mesh, gamma);
  // I.iii Assemble the sparse matrices
  A.setFromTriplets(triplets_A.begin(), triplets_A.end());
  M.setFromTriplets(triplets_M.begin(), triplets_M.end());
  L = A + M;  // Full Galerkin matrix of the LSE \Label[line]{SAApM}
#else
  //====================
  // Your code goes here
  //====================
#endif

  // II. Build the right hand side source vector
#if SOLUTION
  Eigen::VectorXd rhs_vec = computeRHS(mesh, f);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // III. Enforce (zero) dirichlet boundary conditions
#if SOLUTION
  // The B.C. are u(0) = u(N) = 0. The boundary nodes are indexed by 0
  // and N. We can therefore opt for the option of droping first and last rows
  // and columns of L, along with the first and last entries of the rhs\_vec.
  Eigen::SparseMatrix<double> L_reduced = L.block(1, 1, N - 1, N - 1);
  Eigen::VectorXd rhs_vec_reduced = rhs_vec.segment(1, N - 1);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // IV. Solve the LSE L*u = rhs\_vec using an Eigen solver
#if SOLUTION
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  solver.compute(L_reduced);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  u.segment(1, N - 1) = solver.solve(rhs_vec_reduced);
#else
  //====================
  // Your code goes here
  //====================
#endif

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
#if SOLUTION
  std::vector<Eigen::Triplet<double>> triplets = computeA(mesh, alpha);
  A.setFromTriplets(triplets.begin(), triplets.end());
#else
  //====================
  // Your code goes here
  //====================
#endif

  // II. Build the right hand side source vector
#if SOLUTION
  Eigen::VectorXd rhs_vec = computeRHS(mesh, f);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // III. Enforce dirichlet boundary conditions
#if SOLUTION
  // Proceeding as in solveA, we begin by dropping the first and last rows
  // and columns of the galerkin matrix. The rhs\_vec needs to be modified to
  // account for the non-homegenous Dirichlet boundary conditions. It is
  // modified using an offset function.
  Eigen::SparseMatrix<double> A_reduced = A.block(1, 1, N - 1, N - 1);
  Eigen::VectorXd rhs_vec_reduced = rhs_vec.segment(1, N - 1) -
                                    A.block(1, 0, N - 1, 1) * u0 -
                                    A.block(1, N, N - 1, 1) * u1;
#else
  //====================
  // Your code goes here
  //====================
#endif

  // IV. Solve the LSE A*u = rhs\_vec using an built-in Eigen solver
#if SOLUTION
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_reduced);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  u.segment(1, N - 1) = solver.solve(rhs_vec_reduced);
#else
  //====================
  // Your code goes here
  //====================
#endif

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
#if SOLUTION
  // I.i Compute the entries of the Laplace Galerkin matrix A
  std::vector<Eigen::Triplet<double>> triplets_A = computeA(mesh, alpha);
  // I.ii Compute the entries of the mass matrix M
  std::vector<Eigen::Triplet<double>> triplets_M = computeM(mesh, gamma);
  // I.iii Assemble the sparse matrices
  A.setFromTriplets(triplets_A.begin(), triplets_A.end());
  M.setFromTriplets(triplets_M.begin(), triplets_M.end());
  L = A + M;  // Full Galerkin matrix of the LSE
#else
  //====================
  // Your code goes here
  //====================
#endif

  // II. Build the right hand side source vector
#if SOLUTION
  Eigen::VectorXd rhs_vec =
      computeRHS(mesh, [](double) -> double { return 1.0; });
#else
  //====================
  // Your code goes here
  //====================
#endif

  // IV. Solve the LSE A*u = rhs\_vec using an Eigen solver
#if SOLUTION
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(L);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  u = solver.solve(rhs_vec);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return u;
}  // solveC
/* SAM_LISTING_END_C */

}  // namespace LinearFE1D
