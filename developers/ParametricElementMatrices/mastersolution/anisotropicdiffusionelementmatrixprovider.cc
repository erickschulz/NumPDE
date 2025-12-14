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
#if SOLUTION
      // I. OBTAIN COORDINATES OF THE MIDPOINTS
      // I.i Hard-code the midpoints of the edges on the reference triangle
      Eigen::Matrix<double, 2, 3> midpoints_ref(2, 3);
      midpoints_ref << 0.5, 0.5, 0, 0, 0.5, 0.5;
      // I.ii Obtain the midpoints of the parametrized triangle
      auto midpoints_param = cell_geometry->Global(midpoints_ref);

      // II. COMPUTE LOCAL INTEGRATION DATA
      // II.i Compute the inverse of the tranposed Jacobian
      //     (J^T)^{-1} = J*(J^T*J)^{-1}
      // of the parametrization map at the midpoints of the reference triangle
      const Eigen::MatrixXd JinvT(
          cell_geometry->JacobianInverseGramian(midpoints_ref));
      // II.ii Compute the integration element
      //     integration_element(x) = sqrt(det(J^T*J))
      // where J is the Jacobian of the parametrization map
      const Eigen::VectorXd integration_element(
          cell_geometry->IntegrationElement(midpoints_ref));
      // II.iii Hard-code the gradients of the reference basis functions
      Eigen::Matrix<double, 2, 3> gradients_ref(2, 3);
      gradients_ref << -1, 1, 0, -1, 0, 1;

      // III. PERFORM NUMERICAL QUADRATURE
      element_matrix = Eigen::MatrixXd::Zero(3, 3);
      for (int i = 0; i < 3; i++) {  // for each local degree of freedom
        // III.i Evaluate the diffusion tensor
        Eigen::Vector2d anisotropy_vec =
            anisotropy_vec_field_(midpoints_param.col(i));
        Eigen::Matrix2d diffusion_tensor =
            Eigen::Matrix2d::Identity(2, 2) +
            anisotropy_vec * anisotropy_vec.transpose();
        // III. ii Compute gradients of the global basis functions
        Eigen::MatrixXd gradients_param(2, 3);
        gradients_param = JinvT.block(0, 2 * i, 2, 2) * gradients_ref;
        // III. iii Compute local contribution to the element matrix
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            double integrand = (gradients_param.col(j).transpose() *
                                diffusion_tensor * gradients_param.col(k));
            element_matrix(j, k) += integration_element(i) * integrand;
          }
        }
      }
      element_matrix *= 1. / 6.;
#else

      // ===================
      // Your code goes here
      // ===================

#endif
      break;
      /* SAM_LISTING_END_1 */
    }

    /* QUADRILATERAL CELL */
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
      // II.i Compute the inverse of the tranposed Jacobian of the
      //     (J^T)^{-1} = J*(J^T*J)^{-1}
      // parametrization map at the midpoints of the reference quadrilateral
      const Eigen::MatrixXd JinvT(
          cell_geometry->JacobianInverseGramian(midpoints_ref));
      // II.ii Compute the integration element
      //     integration_element(x) = sqrt(det(J^T*J))
      // where J is the Jacobian of the parametrization map
      const Eigen::VectorXd integration_element(
          cell_geometry->IntegrationElement(midpoints_ref));
      // II.iii Hard-code a matrix-valued function that returns the gradients of
      // the reference basis functions evaluated at coordinates
      auto gradients_ref =
          [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 4> {
        Eigen::Matrix<double, 2, 4> element_matrix;
        element_matrix(0, 0) = x(1) - 1;
        element_matrix(1, 0) = x(0) - 1;
        element_matrix(0, 1) = 1 - x(1);
        element_matrix(1, 1) = -x(0);
        element_matrix(0, 2) = x(1);
        element_matrix(1, 2) = x(0);
        element_matrix(0, 3) = -x(1);
        element_matrix(1, 3) = 1 - x(0);
        return element_matrix;
      };

      // III. PERFORM NUMERICAL QUADRATURE
      element_matrix = Eigen::MatrixXd::Zero(4, 4);
      for (int i = 0; i < 4; i++) {  // for each local degree of freedom
        // III.i Evaluate the diffusion tensor
        Eigen::Vector2d anisotropy_vec =
            anisotropy_vec_field_(midpoints_param.col(i));
        Eigen::Matrix2d diffusion_tensor =
            Eigen::Matrix2d::Identity(2, 2) +
            anisotropy_vec * anisotropy_vec.transpose();
        // III.ii Compute gradients of the global basis functions
        Eigen::Matrix<double, 2, 4> gradients_param{
            JinvT.block(0, 2 * i, 2, 2) * gradients_ref(midpoints_ref.col(i))};
        // III.iii Compute local contribution to the element matrix
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            double integrand = (gradients_param.col(j).transpose() *
                                diffusion_tensor * gradients_param.col(k));
            element_matrix(j, k) += integration_element(i) * integrand;
          }
        }
      }
      element_matrix *= 1. / 4.;
#else

      // ===================
      // Your code goes here
      // ===================

#endif
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
