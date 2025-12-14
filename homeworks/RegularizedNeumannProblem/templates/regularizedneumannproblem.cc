/**
 * @file regularizedneumannproblem.cc
 * @brief NPDE homework RegularizedNeumannProblem code
 * @author Christian Mitsch, Philippe Peter
 * @date March 2020
 * @copyright Developed at ETH Zurich
 */

#include "regularizedneumannproblem.h"

namespace RegularizedNeumannProblem {

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd assembleVector_c(const lf::assemble::DofHandler& dofh) {
  Eigen::VectorXd c(dofh.NumDofs());
  //====================
  // Your code goes here
  //====================
  return c;
}
/* SAM_LISTING_END_5 */

}  // namespace RegularizedNeumannProblem
