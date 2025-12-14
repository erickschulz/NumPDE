/**
 * @file sdirk.cc
 * @brief NPDE homework SDIRK code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "sdirk.h"

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.h"

namespace SDIRK {

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d SdirkStep(const Eigen::Vector2d& z0, double h, double gamma) {
  Eigen::Vector2d res;
  // Compute one timestep of the SDIRK implicit RK-SSM for the linear ODE
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> SdirkSolve(const Eigen::Vector2d& z0,
                                        unsigned int M, double T,
                                        double gamma) {
  // Solution vector
  std::vector<Eigen::Vector2d> res(M + 1);
  // Solve the ODE with uniform timesteps using the SDIRK method
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double CvgSDIRK() {
  double conv_rate;
  // Study the convergence rate of the method.
  //====================
  // Your code goes here
  //====================
  return conv_rate;
}
/* SAM_LISTING_END_2 */

}  // namespace SDIRK
