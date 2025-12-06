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

  detJ = h * h * Psi(j * h) + h * h * (Psi((j + 1) * h) - Psi(j * h)) * xhat(0);
  detJ = std::abs(detJ);

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

  invJT(0, 0) = h * Psi(j * h) + h * (Psi((j + 1) * h) - Psi(j * h)) * xhat(0);
  invJT(0, 1) = -h * (Psi((j + 1) * h) - Psi(j * h)) * (l + xhat(1));
  invJT(1, 0) = 0.0;
  invJT(1, 1) = h;

  double detJ = integrationElement(n, j, l, Psi, xhat);
  invJT /= detJ;

  return invJT;
}
/* SAM_LISTING_END_2 */

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

  for (int iq = 0; iq < 4; iq++) {
    // iq-th quadrature node from xq
    Eigen::Vector2d xq_iq = xq.col(iq);
    // Inverse Jacobian transposed at iq-th quadrature node
    Eigen::Matrix2d invJT_iq = jacobianInverseTransposed(n, j, l, Psi, xq_iq);
    // determinant at iq-th quadrature node
    double detJ_iq = integrationElement(n, j, l, Psi, xq_iq);

    // Reference element gradients of basis functions at iq-th node
    Eigen::MatrixXd bhat_grad_iq(2, 4);
    bhat_grad_iq = bhats_grad(xq_iq);

    // Introduce g for a more compact notation, see lecture notes (3.7.36)
    Eigen::MatrixXd g_iq(4, 2);
    g_iq.row(0) = bhat_grad_iq.col(0).transpose() * invJT_iq.transpose();
    g_iq.row(1) = bhat_grad_iq.col(1).transpose() * invJT_iq.transpose();
    g_iq.row(2) = bhat_grad_iq.col(2).transpose() * invJT_iq.transpose();
    g_iq.row(3) = bhat_grad_iq.col(3).transpose() * invJT_iq.transpose();

    // Global coordinates of iq-th quadrature node for alpha
    Eigen::Vector2d xq_gl;
    xq_gl << (j + xq_iq(0)) * h, (l + xq_iq(1)) * h * Psi((j + xq_iq(0)) * h);

    // Weighted contributions to element matrix
    A += wq(iq) * alpha(xq_gl) * detJ_iq * g_iq * g_iq.transpose();
  }

  // For Neumann Boundary Conditions: Surface Contributions
  if (l == n - 1) {
    double dl2 = h * h + (l + 1) * (l + 1) * h * h *
                             (Psi((j + 1) * h) - Psi(j * h)) *
                             (Psi((j + 1) * h) - Psi(j * h));
    double dl = std::sqrt(dl2);

    Eigen::Vector4d bhat2 = bhats(xq.col(2));
    Eigen::Vector4d bhat3 = bhats(xq.col(3));

    A(2, 2) += 0.5 * dl * (bhat2(2) * bhat2(2) + bhat3(2) * bhat3(2));
    A(2, 3) += 0.5 * dl * (bhat2(2) * bhat2(3) + bhat3(2) * bhat3(3));
    A(3, 2) += 0.5 * dl * (bhat2(3) * bhat2(2) + bhat3(3) * bhat3(2));
    A(3, 3) += 0.5 * dl * (bhat2(3) * bhat2(3) + bhat3(3) * bhat3(3));
  }

  return A;
}
/* SAM_LISTING_END_3 */

/* Returns the global index of the local shape function local_dof on element
 * K_jl */
/* SAM_LISTING_BEGIN_4 */
int geoThermLocalToGlobal(unsigned int n, unsigned int j, unsigned int l,
                          unsigned int local_dof) {
  // Map local indices of basis functions to global indices
  int global_dof;

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

  // Loop over elements
  for (int j = 0; j < n; j++) {
    for (int l = 0; l < n; l++) {
      // Obtain element matrix
      Eigen::MatrixXd A_jl = geoThermElemMat(n, j, l, alpha, Psi);

      // Assemble global matrix
      for (int iloc = 0; iloc < 4; iloc++) {
        for (int jloc = 0; jloc < 4; jloc++) {
          // Glonal indices from local dofs
          int iglo = geoThermLocalToGlobal(n, j, l, iloc);
          int jglo = geoThermLocalToGlobal(n, j, l, jloc);

          triplets.push_back(Eigen::Triplet(iglo, jglo, A_jl(iloc, jloc)));
        }
      }
    }
  }

  return triplets;
}
/* SAM_LISTING_END_5 */

/* Replace the m-th row of Galerkin matrix A belonging to any node
 * on the Dirichlet Boundary Gamma_D with the m-th unit vector
 */
/* SAM_LISTING_BEGIN_6 */
void geoThermBdElim(unsigned int n, std::vector<Eigen::Triplet<double>>& A) {
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

  // Get the assembled Galerkin matrix A in triplets format
  std::vector<Eigen::Triplet<double>> triplets =
      assembleGeoTherm(n, alpha, Psi);
  // Eliminate Dirichlet boundary conditions in Galerkin matrix A
  geoThermBdElim(n, triplets);

  // Obtain Sparse Galerkin matrix
  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // Set the right-hand side vector to zero
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N_dofs);
  // Set the right-hand side vector to one accounting for the Dirichlet boundary
  for (int j = 0; j < n + 1; j++) {
    b(j) = 1.0;
  }

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  mu = solver.solve(b);

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

  for (int j = 0; j < n; j++) {
    // Compute length of edfe
    double dl2 = h * h + (Psi((j + 1) * h) - Psi(j * h)) *
                             (Psi((j + 1) * h) - Psi(j * h));
    double dl = std::sqrt(dl2);

    // Trapezoidal rule
    int i = j + (n + 1) * n;
    int ip = j + 1 + (n + 1) * n;
    SurfInt += 0.5 * dl * (mu(i) + mu(ip));
  }

  return SurfInt;
}
/* SAM_LISTING_END_8 */

} /* namespace ParametricFiniteElements */
#endif
