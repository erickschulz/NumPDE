/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "transformationofgalerkinmatrices.h"

#include <cassert>

namespace TransformationOfGalerkinMatrices {

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Triplet<double>> transformCOOmatrix(
    const std::vector<Eigen::Triplet<double>>& A) {
  std::vector<Eigen::Triplet<double>> A_t{};  // return value

  // First step: find the size of the matrix by searching the maximal
  // indices. Depends on the assumption that no zero rows/columns occur.
  int rows_max_idx = 0, cols_max_idx = 0;
  for (const Eigen::Triplet<double>& triplet : A) {
    rows_max_idx =
        (triplet.row() > rows_max_idx) ? triplet.row() : rows_max_idx;
    cols_max_idx =
        (triplet.col() > cols_max_idx) ? triplet.col() : cols_max_idx;
  }
  int n_cols = cols_max_idx + 1;

  // Make sure we deal with a square matrix
  assert(rows_max_idx + 1 == n_cols);
  // The matrix size must have even parity
  assert(n_cols % 2 == 0);

  int N = n_cols;  // Size of (square) matrix
  int M = N / 2;   // Half the size
  //====================
  // Your code goes here
  //====================
  return A_t;
}
/* SAM_LISTING_END_1 */

}  // namespace TransformationOfGalerkinMatrices
