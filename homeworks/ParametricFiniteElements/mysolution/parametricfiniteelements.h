#ifndef PARAMETRICFINITEELEMENTS_HPP
#define PARAMETRICFINITEELEMENTS_HPP
/**
 * @file parametricfiniteelements.h
 * @brief NPDE homework ParametricFiniteElements code
 * @author Amélie Loher
 * @date 04.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <complex>
#include <vector>

namespace ParametricFiniteElements {

/* Computes the absolute value of the determinant of the Jacobian
 * n: mesh discretization parameter
 * j,l: indices of the considered element
 * Psi: topography function
 * xhat: reference element
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
double integrationElement(unsigned int n, unsigned int j, unsigned int l,
                          FUNCTOR&& Psi, Eigen::Vector2d xhat) {
  // Mesh width
  double h = 1.0 / n;
  double detJ = 0.0;

//====================
// Your code goes here
//====================

  return detJ;
}
/* SAM_LISTING_END_1 */

/* Computes the inverse Jacobian transposed
 * n: mesh discretization parameter
 * j,l: indices of the considered element
 * Psi: topography function
 * xhat: reference element
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::Matrix2d jacobianInverseTransposed(unsigned int n, unsigned int j,
                                          unsigned int l, FUNCTOR&& Psi,
                                          Eigen::Vector2d xhat) {
  // Mesh width
  double h = 1.0 / n;
  // Inverse Jacobian transposed
  Eigen::Matrix2d invJT;

//====================
// Your code goes here
//====================
  return invJT;
}
/* SAM_LISTING_END_2 */

/* Returns the basis functions on the Reference Element at node xhat */
inline Eigen::Vector4d bhats(Eigen::Vector2d xhat) {
  Eigen::Vector4d res;

  res(0) = (1 - xhat(0)) * (1 - xhat(1));
  res(1) = xhat(0) * (1 - xhat(1));
  res(2) = xhat(0) * xhat(1);
  res(3) = (1 - xhat(0)) * xhat(1);

  return res;
}

/* Returns the gradients of the basis functions on Reference Element at node
 * xhat */
inline Eigen::MatrixXd bhats_grad(Eigen::Vector2d xhat) {
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

/* Computes the volume contributions to the element matrix for K_j,l
 * with quadrature rule (5.6.6) on problem sheet
 * n: mesh discretization parameter
 * j,l: indices of the considered element
 * Psi: topography function
 * alpha: material function
 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::MatrixXd geoThermElemMat(unsigned int n, unsigned int j, unsigned int l,
                                FUNCTOR1&& alpha, FUNCTOR2&& Psi) {
  // Mesh width
  double h = 1.0 / n;

  // Quadrature nodes and weights on reference element
  Eigen::MatrixXd xq(2, 4);
  xq << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
  Eigen::Vector4d wq;
  wq << 0.25, 0.25, 0.25, 0.25;

  // Volume contributions to element matrix A
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 4);

//====================
// Your code goes here
//====================

  return A;
}
/* SAM_LISTING_END_3 */

/* Returns the global index of the local shape function local_dof on element
 * K_jl */
/* SAM_LISTING_BEGIN_4 */
inline int geoThermLocalToGlobal(unsigned int n, unsigned int j, unsigned int l,
                          unsigned int local_dof) {
  // Map local indices of basis functions to global indices
  int global_dof;

//====================
// Your code goes here
//====================

  return global_dof;
}
/* SAM_LISTING_END_4 */

/* Computes the Galerkin matrix in triplet format based on Element matrix */
/* SAM_LISTING_BEGIN_5 */
template <typename FUNCTOR1, typename FUNCTOR2>
std::vector<Eigen::Triplet<double>> assembleGeoTherm(unsigned int n,
                                                     FUNCTOR1&& alpha,
                                                     FUNCTOR2&& Psi) {
  // Reserve triplets for Galerkin Matrix A
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(4 * 4 * n * n);

//====================
// Your code goes here
//====================

  return triplets;
}
/* SAM_LISTING_END_5 */

/* Replace the m-th row of Galerkin matrix A belonging to any node
 * on the Dirichlet Boundary Gamma_D with the m-th unit vector
 */
/* SAM_LISTING_BEGIN_6 */
inline void geoThermBdElim(unsigned int n, std::vector<Eigen::Triplet<double>>& A) {
//====================
// Your code goes here
//====================
}
/* SAM_LISTING_END_6 */

/* Compute the basis expansion coefficient vector mu of the
 * finite element solution u_N of (5.6.1) in S_1^0(M).
 */
/* SAM_LISTING_BEGIN_7 */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd geoThermSolve(unsigned int n, FUNCTOR1&& alpha,
                              FUNCTOR2&& Psi) {
  // Total Number of dofs
  int N_dofs = (n + 1) * (n + 1);

  // Basis expansion coefficient vector mu
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(N_dofs);

//====================
// Your code goes here
//====================

  return mu;
}
/* SAM_LISTING_END_7 */

/* Computes approximation of the surface integral of u_N over Gamma_S
 * mu: basis expansion coefficients
 * Psi: topography function
 */
/* SAM_LISTING_BEGIN_8 */
template <typename FUNCTOR>
double geoThermSurfInt(unsigned int n, FUNCTOR&& Psi,
                       const Eigen::VectorXd& mu) {
  // Mesh width
  double h = 1.0 / n;

  // Surface integral of geothermal problem
  double SurfInt = 0.0;

//====================
// Your code goes here
//====================

  return SurfInt;
}
/* SAM_LISTING_END_8 */

} /* namespace ParametricFiniteElements */
#endif
