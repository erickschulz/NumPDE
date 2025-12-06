/** @file This homework "TestQuadratureRules" consists of numerically verifying
 * the order of LehrFEM++ quadrature rules
 * @brief NPDE homework TestQuadratureRules
 * @author Erick Schulz, Liaowang Huang (refactoring)
 * @date 08/03/2019, 22/02/2020 (refactoring)
 * @copyright Developed at ETH Zurich
 */

#include <lf/quad/quad.h>

namespace TestQuadratureRules {

/** @brief To do.
 * @param Reference to the QuadRule object to be tested and an unsigned int
 * for the order being tested
 * @return Returns boolean (true if the passed quadrature rule for a triangular
 * reference element has order <order>)
 */
bool testQuadOrderTria(const lf::quad::QuadRule& quad_rule, unsigned int order);
bool testQuadOrderQuad(const lf::quad::QuadRule& quad_rule, unsigned int order);
unsigned int calcQuadOrder(const lf::quad::QuadRule& quad_rule);

}  // namespace TestQuadratureRules
