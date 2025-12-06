#ifndef CLEMPIRICFLUX_H
#define CLEMPIRICFLUX_H

/**
 * @file clempiricflux.h
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

#include "uniformcubicspline.h"

namespace CLEmpiricFlux {

class GodunovFlux {
 public:
  GodunovFlux(const UniformCubicSpline& f);
  // evaluate the Godunov numerical flux F(v, w)
  double operator()(double v, double w) const;

 private:
  // strictly convex flux function (describing the PDE)
  UniformCubicSpline _f;
};

}  // namespace CLEmpiricFlux

#endif
