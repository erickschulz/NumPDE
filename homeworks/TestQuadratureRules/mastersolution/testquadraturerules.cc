/**
 * @file
 * @brief NPDE homework TestQuadratureRules
 * @author Erick Schulz, Liaowang Huang (refactoring)
 * @date 08/03/2019, 22/02/2020 (refactoring)
 * @copyright Developed at ETH Zurich
 */

#include "testquadraturerules.h"

#include <lf/base/base.h>
#include <lf/quad/quad.h>

#include <Eigen/Core>
#include <cassert>
#include <cmath>

namespace TestQuadratureRules {

double factorial(int i) { return std::tgamma(i + 1); }

/* SAM_LISTING_BEGIN_1 */
bool testQuadOrderTria(const lf::quad::QuadRule& quad_rule,
                       unsigned int order) {
  bool order_isExact = true;  // return variable
  double my_epsilon = 1e-12;
  // Retrieve the passed quadrature rule's reference element
  // NOLINTBEGIN(clang-analyzer-deadcode.DeadStores)
  const lf::base::RefEl ref_element = quad_rule.RefEl();
  // NOLINTEND(clang-analyzer-deadcode.DeadStores)
  // Check that the passed reference element is triangular
  assert(ref_element == lf::base::RefElType::kTria);
  // A quadrature rule involves quadrature nodes and weights defined so that
  // the weighted sum of the value of a function at these points approximates
  // the integral of that function.
  const Eigen::VectorXd weights = quad_rule.Weights();
  const Eigen::MatrixXd points = quad_rule.Points();  // (x,y) points
  const Eigen::VectorXd x_coords = points.row(0);
  const Eigen::VectorXd y_coords = points.row(1);
  // A quadrature rule over a two dimensional domain is of order k if it can
  // integrate exactly all bivariate polynomials of order k-1. The collection of
  // such polynomials is spanned by the set of homogeneous polynomials of order
  // strictly less than k, i.e. by the polynomials of the form
  /* p_IJ(x,y) = (x^I)(y^J), I+J < k: */
  auto eval_p_IJ = [&x_coords, &y_coords](int I, int J) -> Eigen::VectorXd {
    return x_coords.array().pow(I) * y_coords.array().pow(J);
  }; /* evaluates p_IJ at all points (x,y) as defined by quad_rule */

  /* Compare analytical value and quadrature sum for all p_IJ, I+J < k */
  double exact_integral;  // analytical value of the integral
  double quad_rule_sum;   // weighted sum used for approximating the integral
  for (int I = 0; I < order; I++) {
    for (int J = 0; J < order - I; J++) {
      exact_integral = factorial(I) * factorial(J) / factorial(I + J + 2);
      quad_rule_sum = eval_p_IJ(I, J).dot(weights);

      // Check if the difference bewteen the results is within tolerance
      order_isExact = fabs(exact_integral - quad_rule_sum) <=
                      fabs(exact_integral) * my_epsilon;
      if (!order_isExact) {
        return order_isExact;
      }
    }
  }
  return order_isExact;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
bool testQuadOrderQuad(const lf::quad::QuadRule& quad_rule,
                       unsigned int order) {
  bool order_isExact = true;  // return variable

  double my_epsilon = 1e-12;
  // Retrieve the passed quadrature rule's reference element
  // NOLINTBEGIN(clang-analyzer-deadcode.DeadStores)
  const lf::base::RefEl ref_element = quad_rule.RefEl();
  // NOLINTEND(clang-analyzer-deadcode.DeadStores)
  // Check that the passed reference element is triangular
  assert(ref_element == lf::base::RefElType::kQuad);
  // A quadrature rule consists of quadrature nodes and weights defined so that
  // the weighted sum of the value of a function at these points approximates
  // the integral of that function.
  const Eigen::VectorXd weights = quad_rule.Weights();
  const Eigen::MatrixXd points = quad_rule.Points();  // (x,y) points
  const Eigen::VectorXd x_coords = points.row(0);
  const Eigen::VectorXd y_coords = points.row(1);
  // A quadrature rule over a two dimensional domain is of order k if it can
  // integrate exactly all bivariate polynomials of order k-1. The collection of
  // such polynomials is spanned by the set of homogeneous polynomials of order
  // strictly less than k, i.e. by the polynomials of the form
  /* p_IJ(x,y) = (1-x)^I(y^J), I,J < k: */
  auto eval_p_IJ = [&x_coords, &y_coords](int I, int J) -> Eigen::VectorXd {
    return x_coords.array().pow(I) * y_coords.array().pow(J);
  }; /* evaluates p_IJ at all points (x,y) as defined by quad_rule */

  /* Compare the analytical value and the quadrature sum for all p_IJ, I+J < k
   */
  double exact_integral;  // analytical value of the integral
  double quad_rule_sum;   // weighted sum used for approximating the integral
  for (int I = 0; I < order; I++) {
    for (int J = 0; J < order; J++) {
      exact_integral = 1.0 / ((I + 1.0) * (J + 1.0));
      quad_rule_sum = eval_p_IJ(I, J).dot(weights);

      // Check if the difference bewteen the results is within tolerance
      order_isExact = fabs(exact_integral - quad_rule_sum) <=
                      fabs(exact_integral) * my_epsilon;
      if (!order_isExact) {
        return order_isExact;
      }
    }
  }
  return order_isExact;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
unsigned int calcQuadOrder(const lf::quad::QuadRule& quad_rule) {
  unsigned int maximal_order = quad_rule.Order();

  // Retrieve the passed quadrature rule's reference element
  const lf::base::RefEl ref_element = quad_rule.RefEl();

  if (ref_element == lf::base::RefElType::kTria) {
    assert(testQuadOrderTria(quad_rule, maximal_order));
    while (testQuadOrderTria(quad_rule, maximal_order + 1)) {
      maximal_order++;
    }
  }

  if (ref_element == lf::base::RefElType::kQuad) {
    assert(testQuadOrderQuad(quad_rule, maximal_order));
    while (testQuadOrderQuad(quad_rule, maximal_order + 1)) {
      maximal_order++;
    }
  }
  return maximal_order;
}
/* SAM_LISTING_END_3 */

}  // namespace TestQuadratureRules
