/**
 * @file clempiricflux.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "clempiricflux.h"

#include <Eigen/Core>
#include <cassert>

namespace CLEmpiricFlux {



GodunovFlux::GodunovFlux(const UniformCubicSpline& f) : _f(f) {};

/* SAM_LISTING_BEGIN_9 */
double GodunovFlux::operator()(double v, double w) const {
  double result;
  //====================
  // Your code goes here
  //====================
  return result;
}

/* SAM_LISTING_END_9 */

}  // namespace CLEmpiricFlux
