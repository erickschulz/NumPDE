/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "evalclass.h"

namespace OutputImpedanceBVP {

/* SAM_LISTING_BEGIN_1 */
EvalResponse::EvalResponse(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p) {
  // Basis vectors for 2D Euclidean space ("unit vectors")
  Eigen::Vector2d e0{1.0, 0.0}, e1{0.0, 1.0};
  // Compute approximate solutions
  Eigen::VectorXd approx_sol_e0, approx_sol_e1;
  approx_sol_e0 = solveImpedanceBVP(fe_space_p, e0);
  approx_sol_e1 = solveImpedanceBVP(fe_space_p, e1);
  // Initialize the
  F_Mat_(0, 0) = computeBoundaryOutputFunctional(approx_sol_e0, fe_space_p, e0);
  F_Mat_(0, 1) = computeBoundaryOutputFunctional(approx_sol_e1, fe_space_p, e0);
  F_Mat_(1, 0) = computeBoundaryOutputFunctional(approx_sol_e0, fe_space_p, e1);
  F_Mat_(1, 1) = computeBoundaryOutputFunctional(approx_sol_e1, fe_space_p, e1);
}

/* Operator overloading */
double EvalResponse::operator()(Eigen::Vector2d g, Eigen::Vector2d d) const {
  double value;
  value = d.dot(F_Mat_ * g);
  return value;
}
/* SAM_LISTING_END_1 */

}  // namespace OutputImpedanceBVP
