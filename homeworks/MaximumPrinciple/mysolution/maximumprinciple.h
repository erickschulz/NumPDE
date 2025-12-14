/**
 * @file  maximumprinciple.h
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>

#ifndef NPDECODES_MAXIMUMPRINCIPLE_MAXIMUMPRINCIPLE_H_
#define NPDECODES_MAXIMUMPRINCIPLE_MAXIMUMPRINCIPLE_H_

namespace MaximumPrinciple {

/**
 * @brief Compute the Galerkin matrix.
 *
 * Computes the global Galerkin Matrix associted with
 * the LHS -(1-c)\Delta u + c u on (0, 1)^2 and the
 * tent functions vanishing at the boundary, for the
 * mesh given in the problem descripition. The
 * local Galerkin matrix of the term c u is computed
 * analytically.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param c Parameter in the equation, with 0<=c<1.
 * @return Galerkin matrix of size M^2 times M^2.
 */
Eigen::SparseMatrix<double> computeGalerkinMatrix(int M, double c);

/**
 * @brief Compute the load vector.
 *
 * Computes the global load vector associted with
 * the source function f on (0, 1)^2 and the tent
 * functions vanishing at the boundary, for the
 * mesh given in the problem descripition.
 * The local integrals are computed by a 2D
 * trapezoidal rule.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param f A functor of type std::function<double(double, double)>.
 * @return Load vector of size M^2.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd computeLoadVector(int M, FUNCTOR&& f) {
  Eigen::VectorXd phi(M * M);
  //====================
  // Your code goes here
  //====================
  return phi;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Compute the Galerkin matrix.
 *
 * Computes the global Galerkin Matrix associted with
 * the LHS -(1-c)\Delta u + c u on (0, 1)^2 and the
 * tent functions vanishing at the boundary, for the
 * mesh given in the problem descripition. The local
 * Galerkin matrix of the term c u is computed by a
 * 2D trapezoidal rule for numerical integration.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param c Parameter in the equation, with 0<=c<1.
 * @return Galerkin matrix of size M^2 times M^2.
 */
Eigen::SparseMatrix<double> computeGalerkinMatrixTR(int M, double c);

}  // namespace MaximumPrinciple

#endif  // NPDECODES_MAXIMUMPRINCIPLE_MAXIMUMPRINCIPLE_H_
