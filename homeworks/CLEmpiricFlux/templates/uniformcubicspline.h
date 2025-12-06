#ifndef UNIFORMCUBICSPLINE_H
#define UNIFORMCUBICSPLINE_H

/**
 * @file uniformcubicspline.h
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>  // Required for assert in helper functions

namespace CLEmpiricFlux {

// Helper functions (moved from uniformcubicspline.cc)
template <typename T>
constexpr T Square(T x) {
  return x * x;
}

template <typename T>
constexpr T Cube(T x) {
  return x * x * x;
}

constexpr int getJ(double a, double b, unsigned int n, double u) {
  return u < b ? (int)(n * ((u - a) / (b - a)) + 1.0) : n;
}

constexpr double zeta(double a, double b, unsigned int n, double j) {
  return a + j * (b - a) / n;
}

class UniformCubicSpline {
 public:
  UniformCubicSpline(double a, double b, Eigen::VectorXd f, Eigen::VectorXd M);
  double operator()(double u) const;  // Point evaluation operator
  double derivative(double u) const;  // Evaluation of derivative
 private:
  unsigned int _n;         // Number of nodes - 1
  double _a, _b;           // Interval boundaries
  Eigen::VectorXd _f;      // Values of flux function at nodes
  Eigen::VectorXd _valsM;  // Values of second derivatives of f at nodes
};

}  // namespace CLEmpiricFlux

#endif
