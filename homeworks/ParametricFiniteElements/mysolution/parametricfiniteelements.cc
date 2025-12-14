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

//====================
// Your code goes here
//====================

  return global_dof;
}
/* SAM_LISTING_END_4 */

/* Replace the m-th row of Galerkin matrix A belonging to any node
 * on the Dirichlet Boundary Gamma_D with the m-th unit vector
 */
/* SAM_LISTING_BEGIN_6 */
void geoThermBdElim(unsigned int n, std::vector<Eigen::Triplet<double>>& A) {
//====================
// Your code goes here
//====================
}
/* SAM_LISTING_END_6 */

}  // namespace ParametricFiniteElements
