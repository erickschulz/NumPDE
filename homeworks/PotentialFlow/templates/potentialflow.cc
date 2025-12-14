#include <Eigen/Core>
#include <Eigen/Sparse>

namespace PotentialFlow {

/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> initializeA(unsigned int M) {
  // For the sake of efficiency the use of Eigen's sparse matrix data type is
  // essential. The matrix is stored in CCS format.
  Eigen::SparseMatrix<double> A(M * M, M * M);
  // We already know that the matrix has at most $9$ non-zero entries per row
  // and column. This information is passed to Eigen via the reserve() member
  // function.
  A.reserve(Eigen::VectorXi::Constant(M * M, 9));
  // Iterate over all interior nodes of the mesh and apply the stencil and
  // initialize the matrix in column-wise order, from top to bottom in every
  // column, which is most efficient for the CCS storage format.
  for (int i = 0; i < M; ++i) {    // "vertical" loop
    for (int j = 0; j < M; ++j) {  // "horizontal" loop
      // Index of the current node
      const int k = i * M + j;
      // Self-interaction weight
      A.insert(k, k) = 16.0 / 6;
      // Interaction term with the node below to the left
      if (i > 0 && j > 0) {
        A.insert(k - M - 1, k) = -2.0 / 6;
      }
      // Interaction term with the node below
      if (i > 0) {
        A.insert(k - M, k) = -2.0 / 6;
      }
      // Interaction term with the node below to the right
      if (i > 0 && j < M - 1) {
        A.insert(k - M + 1, k) = -2.0 / 6;
      }
      // Interaction term with the node to the left
      if (j > 0) {
        A.insert(k - 1, k) = -2.0 / 6;
      }
      // Interaction term with the node to the right
      if (j < M - 1) {
        A.insert(k + 1, k) = -2.0 / 6;
      }
      // Interaction term with the node above to the left
      if (i < M - 1 && j > 0) {
        A.insert(k + M - 1, k) = -2.0 / 6;
      }
      // Interaction term with the node above
      if (i < M - 1) {
        A.insert(k + M, k) = -2.0 / 6;
      }
      // Interaction term with the node above to th right
      if (i < M - 1 && j < M - 1) {
        A.insert(k + M + 1, k) = -2.0 / 6;
      }
    }
  }
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd initializeRHSVector(
    const std::function<double(double, double)>& g, unsigned int M) {
  // Mesh width
  const double h = 1.0 / (M + 1);
  // Off-center entry of stencil
  const double w = 2.0 / 6;
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(M * M);
  // Four corner points
  phi[0] = w * (g(0.0, 0.0) + g(0.0, h) + g(h, 0.0) + g(2 * h, 0.0) +
                g(0.0, 2 * h));  // Bottom left corner
  phi[M - 1] =
      w * (g(1.0, 0.0) + g(1.0 - h, 0.0) + g(1.0, h) + g(1.0 - 2 * h, 0.0) +
           g(1.0, 2 * h));  // Bottom right corner
  phi[M * (M - 1)] =
      w * (g(0.0, 1.0) + g(0.0, 1.0 - h) + g(h, 1.0) + g(2 * h, 1.0) +
           g(0.0, 1.0 - 2 * h));  // Top left corner
  phi[M * M - 1] =
      w * (g(1.0, 1.0) + g(1.0, 1.0 - h) + g(1.0 - h, 1.0) +
           g(1.0, 1.0 - 2 * h) + g(1.0 - 2 * h, 1.0));  // Top right corner
  for (unsigned int l = 1; l < M - 1; ++l) {
    phi[l] = w * (g(h * l, 0.0) + g((l + 1) * h, 0.0) +
                  g((l + 2) * h, 0.0));  // bottom
    phi[M * l] = w * (g(0.0, l * h) + g(0.0, (l + 1) * h) +
                      g(0.0, (l + 2) * h));  // left
    phi[M * (l + 1) - 1] = w * (g(1.0, l * h) + g(1.0, (l + 1) * h) +
                                g(1.0, (l + 2) * h));  // right
    phi[M * (M - 1) + l] =
        w * (g(h * l, 1.0) + g((l + 1) * h, 1.0) + g((l + 2) * h, 1.0));  // top
  }
  return phi;
}
/* SAM_LISTING_END_2 */

}  // namespace PotentialFlow
