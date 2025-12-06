#include "nlmatode.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <iostream>
#include <limits>

// Supplies class Ode45
#include "ode45.h"
// Supplied auxiliary function for linear regression
#include "polyfit.h"

namespace NLMatODE {

/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd& Y0, double T) {
  // Use the Ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  Eigen::MatrixXd YT;
  //====================
  // Your code goes here
  //====================
  return YT;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd& M, double T) {
  // Check if $Y'*Y$ is preserved at the time $T$ by matode.
  //====================
  // Your code goes here
  //====================

  return false;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod() {
  // Compute the fitted convergence rate of the Discrete
  // gradient method. Also tabulate the values M and the errors.
  double conv_rate = 0;
  //====================
  // Your code goes here
  //====================
  return conv_rate;
}
/* SAM_LISTING_END_3 */

}  // namespace NLMatODE
