// Demonstration code for course Numerical Methods for Partial Differential
// Equations Author: R. Hiptmair, SAM, ETH Zurich Date: May 2020

#ifndef CONSFORMEVL_H_
#define CONSFORMEVL_H_

#include <Eigen/Dense>

namespace ConsFV {

using namespace Eigen;

// \ref{cpp:fluxdiff} in lecture notes

// arguments: (Finite) state vector \Blue{$\mu$} of cell averages, see
// \eqref{eq:cellavg}
//  Functor \Blue{$F : \mathbb R \times \mathbb R \mapsto \mathbb R$}, 2-point
//  numerical flux
// return value:  Vector with differences of numerical fluxes, which provides
// the right hand side of \eqref{eq:2pcf}
template <typename FunctionF>
VectorXd fluxdiff(const VectorXd& mu, FunctionF&& F) {
  unsigned n = mu.size();           // length of state vector
  VectorXd fd = VectorXd::Zero(n);  // return vector

  // constant continuation of data for \Blue{$x\leq a$}!
  fd[0] = F(mu[0], mu[1]) - F(mu[0], mu[0]);
  for (unsigned j = 1; j < n - 1; ++j) {
    fd[j] = F(mu[j], mu[j + 1]) - F(mu[j - 1], mu[j]);  // see \eqref{eq:2pcf}
  }
  // constant continuation of data for \Blue{$x\geq b$}!
  fd[n - 1] = F(mu[n - 1], mu[n - 1]) - F(mu[n - 2], mu[n - 1]);
  // Efficient thanks to return value optimization (RVO)
  return fd;
}

}  // namespace ConsFV

#endif
