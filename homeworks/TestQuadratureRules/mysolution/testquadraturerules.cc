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
  //====================
  // Your code goes here
  //====================
  return order_isExact;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
bool testQuadOrderQuad(const lf::quad::QuadRule& quad_rule,
                       unsigned int order) {
  bool order_isExact = true;  // return variable

  //====================
  // Your code goes here
  //====================
  return order_isExact;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
unsigned int calcQuadOrder(const lf::quad::QuadRule& quad_rule) {
  unsigned int maximal_order = quad_rule.Order();

  //====================
  // Your code goes here
  //====================
  return maximal_order;
}
/* SAM_LISTING_END_3 */

}  // namespace TestQuadratureRules
