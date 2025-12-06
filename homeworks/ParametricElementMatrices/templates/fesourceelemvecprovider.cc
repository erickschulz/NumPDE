/** @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans, Erick Schulz (refactoring)
 * @date 13/03/2019, 19/11/2019 (refactoring)
 * @copyright Developed at ETH Zurich */

#include "fesourceelemvecprovider.h"

namespace ParametricElementMatrices {

/** @brief Compute the element vector for
 *
 *       \int_{\Omega} v(x)/(1 + w(x)^2) dx
 *
 * using linear first-order lagrangian finite elements. A local edge-midpoint
 * quadrature rule is used for integration over a cell:
 *
 *   \int_K phi(x) dx = (vol(K)/#vertices(K)) * sum_{edges} phi(midpoint)
 *
 * where K is a cell.
 * @param cell current cell */
Eigen::VectorXd FESourceElemVecProvider::Eval(const lf::mesh::Entity& cell) {
  Eigen::VectorXd element_vector;  // local vector to return;

  /* TOOLS AND DATA */
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_->LocGlobMap()};
  // Obtain cell data
  auto cell_geometry = cell.Geometry();
  auto cell_global_idx = dofh.GlobalDofIndices(cell);

  /* SOLUTION_BEGIN */
  // Integration formula distinguishes between triagular and quadrilateral cells
  switch (cell_geometry->RefEl()) {
    case lf::base::RefEl::kTria(): {
      /* SAM_LISTING_BEGIN_1 */

      // ===================
      // Your code goes here
      // ===================

      break;
      /* SAM_LISTING_END_1 */
    }
    case lf::base::RefEl::kQuad(): {
      /* SAM_LISTING_BEGIN_2 */

      // ===================
      // Your code goes here
      // ===================

      /* SAM_LISTING_END_2 */
    }
      /* ERROR CASE WHERE THE CELL IS NEITHER A TRIANGLE NOR A QUADRILATERAL */
    default:
      LF_VERIFY_MSG(false, "received neither triangle nor quadrilateral");
      /* SOLUTION_END */
  }
  return element_vector;
}
}  // namespace ParametricElementMatrices
