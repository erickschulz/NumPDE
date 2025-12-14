#ifndef SOLVECAUCHYPROBLEM_H
#define SOLVECAUCHYPROBLEM_H

/**
 * @file solvecauchyproblem.h
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 19.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <algorithm>  // For std::min, std::max
#include <cmath>      // For std::abs, std::floor, std::ceil
#include <utility>

#include "clempiricflux.h"
#include "uniformcubicspline.h"

namespace CLEmpiricFlux {

/**
 * @brief Compute an interval containing the support of the solution u at time
 * t.
 *
 * @param f flux function (describing the PDE)
 * @param initsupp interval containing the support of u at initial time
 * @return interval containg the support of u at time t
 */
Eigen::Vector2d findSupport(const UniformCubicSpline& f,
                            Eigen::Vector2d initsupp, double t);

/**
 * @brief Computes the cell averages at initial time on an interval
 * containing the support of the solution u(x,t) for all 0 < t < T.
 *
 * @param f flux function (describing the PDE)
 * @param u0 initial data, -1 < u0 < 1, supported in [-1, 1],
 *                           modeling std::function<double(double)>
 * @param h spacial meshwidth, h > 0.0
 * @param T final time, T > 0.0
 * @return cell averages of u(x, 0)
 */
/* SAM_LISTING_BEGIN_6 */
template <typename FUNCTOR>
Eigen::VectorXd computeInitVec(const UniformCubicSpline& f, FUNCTOR&& u0,
                               double h, double T) {
  Eigen::VectorXd mu0;
  const double A = -1.0;  // left bound of support
  const double B = 1.0;   // right bound of support
  // Here one could also call findSupport()
  const double AT = std::min(A, A + f.derivative(-1.0) * T);
  const double BT = std::max(B, B + f.derivative(1.0) * T);

  const int m_minus = (int)std::floor(AT / h);
  const int m_plus = (int)std::ceil(BT / h);
  const int N = m_plus - m_minus + 1;

  const Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, m_minus * h, m_plus * h);
  mu0 = x.unaryExpr(std::forward<FUNCTOR>(u0));
  return mu0;
}
/* SAM_LISTING_END_6 */

/**
 * @brief Implements the right-hand side of the semi-discretized equation.
 *
 * @param mu0 vector of size N containing the cell averages at initial time
 * @param h spacial mesh-width
 * @param numFlux numerical flux F(v, w) convertible to
 *                                   std::function<double(double, double)>
 * @return vector of size N containg the image of mu0 under the RHS of the
 * semi-discretized equation
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd semiDiscreteRhs(const Eigen::VectorXd& mu0, double h,
                                FUNCTOR&& numFlux) {
  int m = mu0.size();
  Eigen::VectorXd mu1(m);
  mu1(0) = -1.0 / h * (numFlux(mu0(0), mu0(1)) - numFlux(mu0(0), mu0(0)));
  for (int j = 1; j < m - 1; ++j) {
    mu1(j) =
        -1.0 / h * (numFlux(mu0(j), mu0(j + 1)) - numFlux(mu0(j - 1), mu0(j)));
  }
  mu1(m - 1) =
      -1.0 / h *
      (numFlux(mu0(m - 1), mu0(m - 1)) - numFlux(mu0(m - 2), mu0(m - 1)));
  return mu1;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Implements Ralston's method to solve a homogenous ODE
 *
 * @param rhs right-hand side of the homogenous ODE, models
 *                            std::function<Eigen::VectorXd(Eigen::VectorXd)>
 * @param mu0 initial data, vector of size N
 * @param tau timestep size, tau > 0.0
 * @param n number of timesteps to perform, n > 0
 * @return vector of size N containg the approximate solution at time n * tau
 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::VectorXd RalstonODESolver(FUNCTOR&& rhs, Eigen::VectorXd mu0, double tau,
                                 int n) {
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd k1 = rhs(mu0);
    Eigen::VectorXd k2 = rhs(mu0 + tau * 2.0 / 3.0 * k1);
    mu0 = mu0 + 0.25 * tau * (k1 + 3.0 * k2);
  }
  return mu0;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Implements a finite volume scheme to solve a conservation law with
 * strictly convex flux function
 *
 * @param f flux function (describing the PDE)
 * @param mu0 cell averages of initial data
 * @param N number of spacial nodes, N > 1
 * @param T final time, T > 0.0
 * @return cell averages at final time, i.e.of the solution u(x,T)
 */
Eigen::VectorXd solveCauchyProblem(const UniformCubicSpline& f,
                                   const Eigen::VectorXd& mu0, double h,
                                   double T);

}  // namespace CLEmpiricFlux

#endif
