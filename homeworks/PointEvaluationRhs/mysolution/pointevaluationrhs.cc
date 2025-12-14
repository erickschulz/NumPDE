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

  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
  return result;
}

Eigen::Vector2d GlobalInverseQuad(Eigen::Matrix<double, 2, 4> vert,
                                  Eigen::Vector2d x) {
  constexpr double kEPS = 1.0E-8;
  Eigen::Vector2d x_hat;

  // Implement and use the functions triaArea and solveQuadraticEquation
  //====================
  // Your code goes here
  //====================
  return x_hat;
}

std::pair<double, double> normsSolutionPointLoadDirichletBVP(
    const lf::assemble::DofHandler& dofh, Eigen::Vector2d source_point,
    Eigen::VectorXd& sol_vec) {
  std::pair<double, double> result(0, 0);
  const unsigned int N_dofs = dofh.NumDofs();
  sol_vec.resize(N_dofs);
  sol_vec.setZero();
  //====================
  // Your code goes here
  //====================
  return result;
}

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd DeltaLocalVectorAssembler::Eval(const lf::mesh::Entity& cell) {
  Eigen::VectorXd result;
  // get the coordinates of the corners of this cell
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  auto vertices = lf::geometry::Corners(*geo_ptr);
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_6 */

}  // namespace PointEvaluationRhs
