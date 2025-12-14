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
  // Do not forget to initialize vector before assembly!
  c.setZero();
  // ELEMENT_VECTOR_BUILDER object
  VecHelper my_vec_provider_c{};
  // Cell (= codim-0 entities)-oriented assembly into c
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider_c, c);
  return c;
}
/* SAM_LISTING_END_5 */

}  // namespace RegularizedNeumannProblem
