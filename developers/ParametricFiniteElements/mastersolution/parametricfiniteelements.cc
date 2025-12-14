/**
 * @file parametricfiniteelements.cc
 * @brief NPDE homework ParametricFiniteElements code
 * @author Amélie Loher
 * @date 04.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "parametricfiniteelements.h"

#include <Eigen/Core>

namespace ParametricFiniteElements {

/* Returns the basis functions on the Reference Element at node xhat */
Eigen::Vector4d bhats(Eigen::Vector2d xhat) {
  Eigen::Vector4d res;

  res(0) = (1 - xhat(0)) * (1 - xhat(1));
  res(1) = xhat(0) * (1 - xhat(1));
  res(2) = xhat(0) * xhat(1);
  res(3) = (1 - xhat(0)) * xhat(1);

  return res;
}

/* Returns the gradients of the basis functions on Reference Element at node
 * xhat */
Eigen::MatrixXd bhats_grad(Eigen::Vector2d xhat) {
  Eigen::MatrixXd res(2, 4);

  res(0, 0) = xhat(1) - 1;
  res(1, 0) = xhat(0) - 1;
  res(0, 1) = 1 - xhat(1);
  res(1, 1) = -xhat(0);
  res(0, 2) = xhat(1);
  res(1, 2) = xhat(0);
  res(0, 3) = -xhat(1);
  res(1, 3) = 1 - xhat(0);

  return res;
}

/* Returns the global index of the local shape function local_dof on element
 * K_jl */
/* SAM_LISTING_BEGIN_4 */
int geoThermLocalToGlobal(unsigned int n, unsigned int j, unsigned int l,
                          unsigned int local_dof) {
  // Map local indices of basis functions to global indices
  int global_dof;

#if SOLUTION
  switch (local_dof) {
    case 0:
      global_dof = j + (n + 1) * l;
      break;

    case 1:
      global_dof = j + 1 + (n + 1) * l;
      break;

    case 2:
      global_dof = j + 1 + (n + 1) * (l + 1);
      break;

    case 3:
      global_dof = j + (n + 1) * (l + 1);
      break;

    default:
      global_dof = 66;
      break;
  }
#else
//====================
// Your code goes here
//====================
#endif

  return global_dof;
}
/* SAM_LISTING_END_4 */

/* Replace the m-th row of Galerkin matrix A belonging to any node
 * on the Dirichlet Boundary Gamma_D with the m-th unit vector
 */
/* SAM_LISTING_BEGIN_6 */
void geoThermBdElim(unsigned int n, std::vector<Eigen::Triplet<double>>& A) {
#if SOLUTION
  // Identify Triplets on Boundary with Dirichlet Condition
  for (auto& a : A) {
    if (a.row() < n + 1) {
      a = Eigen::Triplet(a.row(), a.col(), 0.0);
    }
  }

  // Set to identity on Dirchlet Boundary part of the boundary Gamma
  for (int i = 0; i < n + 1; i++) {
    A.push_back(Eigen::Triplet(i, i, 1.0));
  }
#else
//====================
// Your code goes here
//====================
#endif
}
/* SAM_LISTING_END_6 */

}  // namespace ParametricFiniteElements
