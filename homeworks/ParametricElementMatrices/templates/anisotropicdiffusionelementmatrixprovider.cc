/** @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans, Erick Schulz (refactoring)
 * @date 13/03/2019, 19/11/2019 (refactoring)
 * @copyright Developed at ETH Zurich */

#include "anisotropicdiffusionelementmatrixprovider.h"

namespace ParametricElementMatrices {

/** @brief Compute the local element matrix for the Galerkin matrix of
 *
 *     \int_{\Omega} (1 + d(x)d(x)') * grad(u(x)).grad(v(x)) dx
 *
 * using linear first-order lagrangian finite elements. A local edge-midpoint
 * quadrature rule is used for integration over a cell:
 *
 *   \int_K phi(x) dx = (vol(K)/#vertices(K)) * sum_{edges} phi(midpoint)
 *
 * where K is a cell.
 * @param cell current cell */
Eigen::MatrixXd AnisotropicDiffusionElementMatrixProvider::Eval(
    const lf::mesh::Entity& cell) {
  Eigen::MatrixXd element_matrix;  // local matrix to return

  // Cell data
  auto cell_geometry = cell.Geometry();

  /* SOLUTION_BEGIN */
  // Integration formula distinguishes between triagular and quadrilateral cells
  switch (cell_geometry->RefEl()) {
    /* TRIANGULAR CELL */
    case lf::base::RefEl::kTria(): {
      /* SAM_LISTING_BEGIN_1 */

      // ===================
      // Your code goes here
      // ===================

      break;
      /* SAM_LISTING_END_1 */
    }

    /* QUADRILATERAL CELL */
    case lf::base::RefEl::kQuad(): {
      /* SAM_LISTING_BEGIN_2 */

      // ===================
      // Your code goes here
      // ===================

      break;
      /* SAM_LISTING_END_2 */
    }

    /* ERROR CASE WHERE THE CELL IS NEITHER A TRIANGLE NOR A QUADRILATERAL */
    default:
      LF_VERIFY_MSG(false, "received neither triangle nor quadrilateral");
  }
  /* SOLUTION_END */
  return element_matrix;
}  // AnisotropicDiffusionElementMatrixProvider::Eval

}  // namespace ParametricElementMatrices
