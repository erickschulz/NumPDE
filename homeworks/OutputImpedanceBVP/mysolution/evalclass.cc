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
  //====================
  // Your code goes here
  //====================
}

/* Operator overloading */
double EvalResponse::operator()(Eigen::Vector2d g, Eigen::Vector2d d) const {
  double value;
  //====================
  // Your code goes here
  //====================
  return value;
}
/* SAM_LISTING_END_1 */

}  // namespace OutputImpedanceBVP
