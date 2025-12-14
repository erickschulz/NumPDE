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
  if (v >= w)
    result = std::max(_f(v), _f(w));
  else {
    if (_f.derivative(v) > 0.0)
      result = _f(v);
    else if (_f.derivative(w) < 0.0)
      result = _f(w);
    else {
      auto df = [this](double x) { return _f.derivative(x); };
      double z = findRoots(v, w, df);
      result = _f(z);
    }
  }
  return result;
}

/* SAM_LISTING_END_9 */

}  // namespace CLEmpiricFlux
