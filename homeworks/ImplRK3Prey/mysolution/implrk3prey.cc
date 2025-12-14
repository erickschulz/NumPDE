/**
 * @file implrk3prey.cc
 * @brief NPDE homework ImplRK3Prey code
 * @author Unknown, Oliver Rietmann
 * @date 29.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "implrk3prey.h"

#include <Eigen/Core>

namespace ImplRK3Prey {

// Compute the Kronecker product $C = A \otimes B$
// A is m x n matrix, B is l x k matrix
// return Kronecker product of A and B: dim is m*l x n*k
Eigen::MatrixXd kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }
  return C;
}

}  // namespace ImplRK3Prey
