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
#if SOLUTION
  double near_neighbour_contribution = 2.0 * B_K(0, 1);
  double far_neighbour_contribution = 2.0 * B_K(1, 2);
  double self_contribution = 2.0 * (B_K(0, 0) + B_K(1, 1) + B_K(2, 2));

  std::vector<double> contribution = {
      far_neighbour_contribution,  near_neighbour_contribution,
      near_neighbour_contribution, self_contribution,
      near_neighbour_contribution, near_neighbour_contribution,
      far_neighbour_contribution};

  std::vector<Eigen::Vector2i> shift = {{-1, -1}, {0, -1}, {-1, 0}, {0, 0},
                                        {1, 0},   {0, 1},  {1, 1}};

  std::vector<Eigen::Triplet<double>> tripletList;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      Eigen::Vector2i self = Eigen::Vector2i(i, j);
      for (int k = 0; k < 7; ++k) {
        Eigen::Vector2i other = self + shift[k];
        if (0 <= other(0) && other(0) < M && 0 <= other(1) && other(1) < M) {
          tripletList.push_back(Eigen::Triplet<double>(
              self(0) + M * self(1), other(0) + M * other(1), contribution[k]));
        }
      }
    }
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
#else
  //====================
  // Your code goes here
  //====================
#endif
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double> computeGalerkinMatrix(int M, double c) {
  Eigen::Matrix3d B_K;
#if SOLUTION
  double h = 1.0 / (M + 1);
  Eigen::Matrix3d A_K;
  A_K << 1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5;
  Eigen::Matrix3d M_K;
  M_K << 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0;
  M_K *= h * h / 24.0;
  B_K = (1.0 - c) * A_K + c * M_K;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return assemble(M, B_K);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<double> computeGalerkinMatrixTR(int M, double c) {
  Eigen::Matrix3d B_K;
#if SOLUTION
  double h = 1.0 / (M + 1);
  Eigen::Matrix3d A_K;
  A_K << 1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5;
  Eigen::Matrix3d M_K = h * h / 6.0 * Eigen::Matrix3d::Identity();
  B_K = (1.0 - c) * A_K + c * M_K;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return assemble(M, B_K);
}
/* SAM_LISTING_END_4 */

}  // namespace MaximumPrinciple
