/**
 * @ file pointEvaluation.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch, Liaowang Huang (refactoring)
 * @ date 22/03/2019, 06/01/2020 (refactoring)
 * @ copyright Developed at ETH Zurich
 */

#include "pointevaluationrhs.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <cmath>
#include <utility>

#include "pointevaluationrhs_norms.h"

namespace PointEvaluationRhs {

/* SAM_LISTING_BEGIN_5 */
Eigen::Vector2d GlobalInverseTria(Eigen::Matrix<double, 2, 3> mycorners,
                                  Eigen::Vector2d x) {
  Eigen::Vector2d x_hat;

#if SOLUTION
  // The unique affine mapping of the unit triangle to a general triangle
  // is given by the formula $\cob{\Bx = \VA_K*\wh{\Bx} + \Bs}$,
  // Use the vertex coordinates to calculate matrix A and translation vector s
  Eigen::Matrix2d A(2, 2);
  A.col(0) = mycorners.col(1) - mycorners.col(0);
  A.col(1) = mycorners.col(2) - mycorners.col(0);
  // The inverse mapping: $\cob{\wh{\Bx} = \VA_K^{-1} * (\Bx - \Ba_0)}$
  x_hat = A.partialPivLu().solve(x - mycorners.col(0));
#else
  //====================
  // Your code goes here
  //====================
#endif
  return x_hat;
}
/* SAM_LISTING_END_5 */

/** @brief Numerically stable solution of a quadratic equation in R
 * @param a,b,c coefficients of quadratic polynomias ax^2+bx+c
 * @return both zeros, NaN if complex
 */
/* SAM_LISTING_BEGIN_4 */
std::pair<double, double> solveQuadraticEquation(double a, double b, double c) {
  // Implement the cases which are solvable and return their solutions
#if SOLUTION
  if (a != 0.) {
    b /= a;
    c /= a;
    const double D = b * b - 4 * c;  // discriminant
    if (D >= 0) {                    // real solutions
      const double rtD = std::sqrt(D);
      // Real solutions, cancellation-free formulas !
      if (b < 0) {
        const double root = 0.5 * (-b + rtD);
        return {root, c / root};
      }  // b >= 0
      const double root = 0.5 * (-b - rtD);
      return {root, c / root};
    }
  } else {
    if (b != 0.0) {
      return {-c / b, -c / b};
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  // Return NAN if there are no (real) roots
  return {NAN, NAN};
}

/* SAM_LISTING_END_4 */

/** @brief Computes the area of a triangle
 * @param a,b,c vertex coordinate vectors
 */
inline double triaArea(const Eigen::Vector2d a, const Eigen::Vector2d b,
                       const Eigen::Vector2d c) {
  double result = 0;
#if SOLUTION
  result = 0.5 * std::abs((b[0] - a[0]) * (c[1] - a[1]) -
                          (b[1] - a[1]) * (c[0] - a[0]));
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}

Eigen::Vector2d GlobalInverseQuad(Eigen::Matrix<double, 2, 4> vert,
                                  Eigen::Vector2d x) {
  constexpr double kEPS = 1.0E-8;
  Eigen::Vector2d x_hat;

  // Implement and use the functions triaArea and solveQuadraticEquation
#if SOLUTION
  /* SAM_LISTING_BEGIN_1 */
  // I: Find sub-triangle with largest area and use its
  // middle vertex as new "vertex 0".
  unsigned int vt_zero_idx{0};
  double max_area{0.0};
  for (int l = 0; l < 4; ++l) {
    const double sub_tria_area =
        triaArea(vert.col(l), vert.col((l + 1) % 4), vert.col((l + 2) % 4));
    if (sub_tria_area > max_area) {
      vt_zero_idx = (l + 1) % 4;
      max_area = sub_tria_area;
    }
  }
  // vt_zero_idx stores the index of the vertex that will now
  // be regarded as "pivotal vertex 0"

  /* SAM_LISTING_END_1 */

  const Eigen::Vector2d corner0 = vert.col(vt_zero_idx);
  const Eigen::Vector2d corner1 = vert.col((vt_zero_idx + 1) % 4);
  const Eigen::Vector2d corner2 = vert.col((vt_zero_idx + 2) % 4);
  const Eigen::Vector2d corner3 = vert.col((vt_zero_idx + 3) % 4);

  /* SAM_LISTING_BEGIN_2 */
  // Use corners to calculate matrix A, the translation t
  // and the coefficient vector d
  Eigen::Matrix2d A(2, 2);
  A.col(0) = corner1 - corner0;
  A.col(1) = corner3 - corner0;
  Eigen::Vector2d t = corner0;
  Eigen::Vector2d d = corner2 + corner0 - corner1 - corner3;

  // Calculate the coefficients in the non-linear system of equations
  Eigen::Vector2d y = A.partialPivLu().solve(x - t);
  Eigen::Vector2d q = A.partialPivLu().solve(d);
  Eigen::Vector2d q_orth(q(1), -q(0));

  // Treat exceptional cases
  // If q vanishes
  const double ref_size = y.norm();
  if (q.norm() < kEPS * ref_size) {
    x_hat = y;
  } else if (std::abs(q(0)) < kEPS * ref_size) {
    x_hat(0) = y(0);
    x_hat(1) = y(1) / (1 + q(1) * x_hat(0));
  } else if (std::abs(q(1)) < kEPS * ref_size) {
    x_hat(1) = y(1);
    x_hat(0) = y(0) / (1 + q(0) * x_hat(1));
  } else {
    // Generic case; solve quadratic equation
    double a = q(0);
    double b = 1 + y.dot(q_orth);
    double c = -y(1);

    auto [x_op1, x_op2] = solveQuadraticEquation(a, b, c);
    // No solution
    if (std::isnan(x_op1) || std::isnan(x_op2)) {
      x_hat(0) = NAN;
      x_hat(1) = NAN;
    } else {
      // There is a solution
      // Find root in the unit interval -> x_op1
      if ((x_op1 < 0.) || (x_op1 > 1.)) {
        if ((x_op2 >= 0.) && (x_op2 <= 1.)) {
          x_op1 = x_op2;
        }
      }
      x_hat(1) = x_op1;
      x_hat(0) = (y.dot(q_orth) + q(0) * x_hat(1)) / q(1);
    }
  }
  /* SAM_LISTING_END_2 */
  /* SAM_LISTING_BEGIN_3 */
  // Reverse initial permutation
  switch (vt_zero_idx) {
    case 1: {
      x_hat = ((Eigen::Matrix2d() << 0., -1., 1., 0.)).finished() * x_hat +
              Eigen::Vector2d(1.0, 0.0);
      break;
    }
    case 2: {
      x_hat = -x_hat + Eigen::Vector2d(1.0, 1.0);
      break;
    }
    case 3: {
      x_hat = ((Eigen::Matrix2d() << 0., 1., -1., 0.)).finished() * x_hat +
              Eigen::Vector2d(0.0, 1.0);
      break;
    }
  }
    /* SAM_LISTING_END_3 */
#else
  //====================
  // Your code goes here
  //====================
#endif
  return x_hat;
}

std::pair<double, double> normsSolutionPointLoadDirichletBVP(
    const lf::assemble::DofHandler& dofh, Eigen::Vector2d source_point,
    Eigen::VectorXd& sol_vec) {
  std::pair<double, double> result(0, 0);
  const unsigned int N_dofs = dofh.NumDofs();
  sol_vec.resize(N_dofs);
  sol_vec.setZero();
#if SOLUTION
  /* SAM_LISTING_BEGIN_7 */
  // Assemble matrix A
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix loc_mat_laplace{};
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, loc_mat_laplace, A);
  // Build rhs vector using dedicated ENTITY_VECTOR_PROVIDER
  Eigen::VectorXd rhs(N_dofs);
  rhs.setZero();
  PointEvaluationRhs::DeltaLocalVectorAssembler myvec_pro(source_point);
  lf::assemble::AssembleVectorLocally(0, dofh, myvec_pro, rhs);

  // Enforce Dirichlet boundary conditions
  const double boundary_val = 0;  // zero Dirichlet boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2)};
  auto my_selector = [&dofh, &bd_flags, &boundary_val](unsigned int dof_idx) {
    if (bd_flags(dofh.Entity(dof_idx))) {
      return (std::pair<bool, double>(true, boundary_val));
    }  // interior node: the value we return here does not matter
    return (std::pair<bool, double>(false, 42.0));
  };
  lf::assemble::FixFlaggedSolutionComponents<double>(my_selector, A, rhs);

  // Solve linear system of equations A*x = rhs
  const Eigen::SparseMatrix<double> A_crs(A.makeSparse());
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  if (solver.info() == Eigen::Success) {
    sol_vec = solver.solve(rhs);
  } else {
    LF_ASSERT_MSG(false, "Eigen Factorization failed");
  }
  /* SAM_LISTING_END_7 */
  // return the norms of the solution vector
  result = std::pair<double, double>(
      PointEvaluationRhs::computeL2normLinearFE(dofh, sol_vec),
      PointEvaluationRhs::computeH1seminormLinearFE(dofh, sol_vec));
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd DeltaLocalVectorAssembler::Eval(const lf::mesh::Entity& cell) {
  Eigen::VectorXd result;
  // get the coordinates of the corners of this cell
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  auto vertices = lf::geometry::Corners(*geo_ptr);
#if SOLUTION
  Eigen::Vector2d x_hat;
  // the margin we allow when we determine wether a point is inside an
  // element
  const double margin = 1e-10;
  if (lf::base::RefEl::kTria() == cell.RefEl()) {
    x_hat = PointEvaluationRhs::GlobalInverseTria(vertices, x_0);
    result.resize(3);
    result.setZero();
    if (x_hat(0) <= 1 + margin && x_hat(0) >= 0 - margin &&
        x_hat(1) <= 1 + margin && x_hat(1) >= 0 - margin &&
        x_hat(1) + x_hat(0) <= 1 + margin) {
      already_found = true;
      // Barycentric coordinates on reference triangle
      result[0] = 1.0 - x_hat(0) - x_hat(1);
      result[1] = x_hat(0);
      result[2] = x_hat(1);
    }
  } else if (lf::base::RefEl::kQuad() == cell.RefEl()) {
    x_hat = PointEvaluationRhs::GlobalInverseQuad(vertices, x_0);
    result.resize(4);
    result.setZero();
    if (x_hat(0) <= 1 + margin && x_hat(0) >= 0 - margin &&
        x_hat(1) <= 1 + margin && x_hat(1) >= 0 - margin) {
      already_found = true;
      // Local shape functions on unit square
      result[0] = (1.0 - x_hat(0)) * (1.0 - x_hat(1));
      result[1] = x_hat(0) * (1.0 - x_hat(1));
      result[2] = x_hat(0) * x_hat(1);
      result[3] = (1.0 - x_hat(0)) * x_hat(1);
    }
  } else {
    LF_ASSERT_MSG(
        false, "Function only defined for triangular or quadrilateral cells");
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_6 */

}  // namespace PointEvaluationRhs
