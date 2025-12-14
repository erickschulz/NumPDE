/** @file
This homework "TransformationOfGalerkinMatrices" consists of matrices
transformation (change of Galerkin bases) using the triplet format.
 */

#include <Eigen/Sparse>
#include <vector>

namespace TransformationOfGalerkinMatrices {

/** @brief Transformation of Galerking matrix in COO format using the change of
 * basis described in the problem statement.
 * @param A "Old" Galerkin matrix in COO format
 * @return triplets describing "New" Galerkin matrix
 */
std::vector<Eigen::Triplet<double>> transformCOOmatrix(
    const std::vector<Eigen::Triplet<double>>& A);

}  // namespace TransformationOfGalerkinMatrices
