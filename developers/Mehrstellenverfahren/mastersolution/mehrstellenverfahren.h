/**
 * @ file
 * @ brief NPDE homework TEMPLATE HEADER FILE
 * @ author Tobias Rohner
 * @ date 25-03-2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <iostream>

namespace mehrstellenverfahren {

/**
 * @brief Compute the stiffness matrix A using the Mehrstellenverfahren
 * @param M Size of the tensor product mesh
 * @returns Stiffness matrix A
 */
Eigen::SparseMatrix<double> compMehrstellenA(unsigned int M);

/**
 * @brief Compute the load vector phi using the Mehrstellenverfahren
 * @param f Right hand side function f taking two doubles as x- and
 * y-coordinates
 * @param M Size of the tensor product mesh
 * @returns The load vector phi
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
Eigen::VectorXd compMehrstellenf(FUNCTOR&& f, unsigned int M) {
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(M * M);
#if SOLUTION
  const double h = 1. / (M + 1);
  // Iterate over all interior nodes of the mesh
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      // Translate i and j to x- and y-coordinates
      const double x = (i + 1) * h;
      const double y = (j + 1) * h;
      // Compute the index of the node at (i,j)
      const int k = i * M + j;
      // Compute the k-th element of the load vector
      phi[k] += f(x, y - h);
      phi[k] += f(x - h, y);
      phi[k] += 8 * f(x, y);
      phi[k] += f(x + h, y);
      phi[k] += f(x, y + h);
    }
  }
  phi *= h * h / 12;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return phi;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solve the Poisson equation using the Mehrstellenverfahren
 * @param f Right hand side function f taking two doubles as x- and
 * y-coordinates
 * @param M Size of the tensor product mesh
 * @returns The solution to the Poisson problem
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd solveMehrstellen(FUNCTOR&& f, unsigned int M) {
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(M * M);
#if SOLUTION
  // Compute the stiffness matrix
  Eigen::SparseMatrix<double> A = compMehrstellenA(M);
  // Compute the load vector
  Eigen::VectorXd phi = compMehrstellenf(f, M);
  // Solve the LSE using sparse cholesky
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  mu = solver.solve(phi);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return mu;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Compute the error at a given mesh resolution
 * @param M Size of the tensor product mesh
 * @returns The error of the Mehrstellen discretization
 */
double compgriderr(unsigned int M);

/**
 * @brief Prints out a table containing the error for the Mehrstellenverfahren
 * for various mesh resolutions
 */
void tabulateMehrstellenError();

}  // namespace mehrstellenverfahren
