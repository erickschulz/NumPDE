/* **********************************************************************
 * Demonstration code for NPDE lecture & homeworks
 ********************************************************************** */

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace std;

int main(int, char**) {
  cout << "Demonstration of initialization of sparse matrix in eigen" << endl;

  /* SAM_LISTING_BEGIN_1 */
  const int n = 20, m = 10;
  // Set up zero sparse matrix with row major storage format
  // This format is essential for being able to set the maximal
  // number of non-zero entries \textbf{per row}.
  Eigen::SparseMatrix<int, Eigen::RowMajor> X(n, m);
  // Reserve space for at most nnz\_row non-zero entries per row
  const std::size_t nnz_row = 3;
  X.reserve(Eigen::VectorXi::Constant(n, nnz_row));
  // Initialize nnz\_row  entries per row
  for (int row_idx = 0; row_idx < n; ++row_idx) {
    for (int k = 0; k < nnz_row; ++k) {
      X.coeffRef(row_idx, (row_idx * k) % m) += 1;
    }
  }
  /* SAM_LISTING_END_1 */
  Eigen::MatrixXi X_dense = X;
  cout << "Matrix X = " << endl << X_dense << endl;
  return 0;
}
