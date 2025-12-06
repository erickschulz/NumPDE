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
  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
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
