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
#if SOLUTION
      // I. OBTAIN COORDINATES OF THE MIDPOINTS
      // I.i Hard-code the midpoints of the edges on the reference triangle
      Eigen::MatrixXd midpoints_ref(2, 3);
      midpoints_ref << 0.5, 0.5, 0, 0, 0.5, 0.5;
      // I.ii Obtain the midpoints of the parametrized triangle
      auto midpoints_param = cell_geometry->Global(midpoints_ref);

      // II. COMPUTE LOCAL INTEGRATION DATA
      // II.i Compute the integration element
      //     integration_element(x) = sqrt(det(J^T*J))
      // where J is the Jacobian of the parametrization map
      const Eigen::VectorXd integration_element(
          cell_geometry->IntegrationElement(midpoints_ref));
      // II.ii Evaluate w(x) at the midpoints
      Eigen::VectorXd w(3);
      for (int i = 0; i < 3; i++) {
        w(i) = coeff_expansion_(cell_global_idx[i]);
      }

      // III. PERFORM NUMERICAL QUADRATURE
      element_vector = Eigen::VectorXd::Zero(3);
      for (int i = 0; i < 3; i++) {  // for each local degree of freedom
        element_vector(i) +=
            integration_element(i) *
            (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 3), 2)));
        element_vector((i + 1) % 3) +=
            integration_element(i) *
            (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 3), 2)));
      }
      element_vector *= 1. / 6.;
#else

      // ===================
      // Your code goes here
      // ===================

#endif
      break;
      /* SAM_LISTING_END_1 */
    }
    case lf::base::RefEl::kQuad(): {
      /* SAM_LISTING_BEGIN_2 */
#if SOLUTION
      // I. OBTAIN COORDINATES OF THE MIDPOINTS
      // I.i Hard-code the midpoints of the edges on the reference quadrilateral
      Eigen::Matrix<double, 2, 4> midpoints_ref(2, 4);
      midpoints_ref << 0.5, 1, 0.5, 0, 0, 0.5, 1, 0.5;
      // I.ii Obtain the midpoints of the parametrized quadrilateral
      auto midpoints_param = cell_geometry->Global(midpoints_ref);

      // II. COMPUTE LOCAL INTEGRATION DATA
      // II.i Compute the integration element
      //     integration_element(x) = sqrt(det(J^T*J))
      // where J is the Jacobian of the parametrization map
      const Eigen::VectorXd integration_element(
          cell_geometry->IntegrationElement(midpoints_ref));
      // II.ii Evaluate w(x) at the midpoints
      Eigen::VectorXd w(4);
      for (int i = 0; i < 4; i++) {
        w(i) = coeff_expansion_(cell_global_idx[i]);
      }

      // III. PERFORM NUMERICAL QUADRATURE
      element_vector = Eigen::VectorXd::Zero(4);
      for (int i = 0; i < 4; i++) {  // for each local degree of freedom
        element_vector(i) +=
            integration_element(i) *
            (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 4), 2)));
        element_vector((i + 1) % 4) +=
            integration_element(i) *
            (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 4), 2)));
      }
      element_vector *= 1. / 4.;
      break;
#else

      // ===================
      // Your code goes here
      // ===================

#endif
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
