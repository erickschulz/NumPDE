/**
 * @file finitevolumesineconslaw.cc
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "finitevolumesineconslaw.h"

#include <Eigen/Core>
#include <cmath>

namespace FiniteVolumeSineConsLaw {

/* SAM_LISTING_BEGIN_1 */
constexpr double PI = 3.14159265358979323846;

double f(double x) { return std::sin(PI * x); }

double sineGodFlux(double v, double w) {
  double result;
  //====================
  // Your code goes here
  //====================

  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd sineClawRhs(const Eigen::VectorXd& mu) {
  int N = mu.size();
  double h = 12.0 / N;
  Eigen::VectorXd result(N);

  //====================
  // Your code goes here
  //====================

  return result;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
bool blowup(const Eigen::VectorXd& mu) {
  return mu.minCoeff() < 0.0 || mu.maxCoeff() > 2.0;
}

unsigned int findTimesteps() {
  const unsigned int N = 600;
  unsigned int ML = 100;
  unsigned int MR = 200;

  //====================
  // Your code goes here
  //====================

  return MR;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd sineClawReactionRhs(const Eigen::VectorXd& mu, double c) {
  Eigen::VectorXd rhs(mu.size());
  //====================
  // Your code goes here
  //====================
  return rhs;
}
/* SAM_LISTING_END_4 */

}  // namespace FiniteVolumeSineConsLaw
