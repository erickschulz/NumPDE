/**
 * @file  maximumprinciple.cc
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "maximumprinciple.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

namespace MaximumPrinciple {

/**
 * @brief Assembly on a tensor product mesh
 *
 * Compute the global Galerkin matrix from the local
 * element matrix.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param B_K Local element matrix.
 * @return Global Galerkin matrix of size M^2 times M^2.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> assemble(int M, const Eigen::Matrix3d& B_K) {
  int M2 = M * M;
  Eigen::SparseMatrix<double> A(M2, M2);
  //====================
  // Your code goes here
  //====================
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double> computeGalerkinMatrix(int M, double c) {
  Eigen::Matrix3d B_K;
  //====================
  // Your code goes here
  //====================
  return assemble(M, B_K);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<double> computeGalerkinMatrixTR(int M, double c) {
  Eigen::Matrix3d B_K;
  //====================
  // Your code goes here
  //====================
  return assemble(M, B_K);
}
/* SAM_LISTING_END_4 */

}  // namespace MaximumPrinciple
