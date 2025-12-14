#ifndef EVALCLASS_HPP
#define EVALCLASS_HPP

/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "outputimpedancebvp.h"

namespace OutputImpedanceBVP {

/* SAM_LISTING_BEGIN_1 */
class EvalResponse {
 public:
  /* Constructor */
  explicit EvalResponse(
      const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>&
          fe_space_p);
  /* Evaluation operator */
  double operator()(Eigen::Vector2d g, Eigen::Vector2d d) const;

 private:
  //====================
  // Your code goes here
  //====================
};  // class EvalResponse
/* SAM_LISTING_END_1 */

}  // namespace OutputImpedanceBVP

#endif
