/**
 * @file fvpsystem.h
 * @brief NPDE homework FVPsystem code
 * @ author Wouter Tonnon
 * @ date May 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef FVPsystem_H_
#define FVPsystem_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <Eigen/Core>
#include <iostream>

namespace FVPsystem {

/** @brief Fully discrete finite volume method for the p-system based on HLLE
 * numerical flux
 *
 * @tparam u0Functor std::function<Eigen::Vector2d(double)>
 * @tparam recorder std::function<void(const Eigen::MatrixXd &)
 * @param a left bound of interval of interest
 * @param b right bound of interval of interest
 * @param T final time
 * @param N number of mesh cell in interval of interest
 * @param M number of uniform timesteps
 * @param u0 initial data
 *
 * Special case p(v) = -exp(v)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename u0Functor,
          typename RECORDER = std::function<void(const Eigen::MatrixXd&)>>
Eigen::MatrixXd ev1ExpPSystem(
    double a, double b, double T, unsigned int N, unsigned int M,
    u0Functor&& u0, RECORDER recorder = [](const Eigen::MatrixXd&) -> void {}) {
  //  Determine constants of scheme
  const double dt = T / M;       // timestep size
  const double h = (b - a) / N;  // meshwidth
  Eigen::MatrixXd mu(2, N);      // return vector

#if SOLUTION
  // Set initial conditions
  for (unsigned j = 0; j < N; ++j) {
    mu.col(j) = u0(a + h / 2. + j * h);
  }

  // Define the function p, note that $p' = p$ for this special case
  auto p = [](double u) -> double { return -std::exp(u); };

  // Define the continuous flux function $\VF$, see \prbeqref{eq:psysu}
  auto F = [p](Eigen::Vector2d v) -> Eigen::Vector2d {
    return {-v[1], p(v[0])};
  };

  // Define the HLLE flux \lref{eq:HLLEnfs}
  auto numflux = [p, F](Eigen::MatrixXd v,
                        Eigen::MatrixXd w) -> Eigen::VectorXd {
    double smin = std::min(-std::sqrt(-p(v(0))),
                           -std::sqrt(-(p(w(0)) - p(v(0))) / (w(0) - v(0))));
    double smax = std::max(std::sqrt(-p(w(0))),
                           std::sqrt(-(p(w(0)) - p(v(0))) / (w(0) - v(0))));
    if (smin > 0) return F(v);
    if (smax > 0) {
      const Eigen::Vector2d ustar =
          1.0 / (smin - smax) * (F(w) - F(v) - smax * w + smin * v);
      return F(ustar);
    }
    return F(w);
  };

  // Solve the MOL ODE using explicit Euler timestepping
  recorder(mu);
  for (double t = 0; t < T - dt / 2.; t += dt) {
    Eigen::MatrixXd fd(2, mu.cols());  // Flux difference

    // constant continuation of data
    fd.col(0) = numflux(mu.col(0), mu.col(1)) - numflux(mu.col(0), mu.col(0));
    // Fluxes on the interior
    for (unsigned j = 1; j < N - 1; ++j) {
      fd.col(j) = numflux(mu.col(j), mu.col(j + 1)) -
                  numflux(mu.col(j - 1), mu.col(j));  // see \eqref{eq:2pcf}
    }
    // constant continuation of data!
    fd.col(N - 1) = numflux(mu.col(N - 1), mu.col(N - 1)) -
                    numflux(mu.col(N - 2), mu.col(N - 1));

    mu -= dt / h * fd;
    recorder(mu);
  }
#else
  /* Your code goes here! */
#endif

  // Return
  return mu;
}
/* SAM_LISTING_END_1 */

}  // namespace FVPsystem

#endif
