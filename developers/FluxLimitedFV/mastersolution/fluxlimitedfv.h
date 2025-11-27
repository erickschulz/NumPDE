/**
 * @ file
 * @ brief NPDE homework TEMPLATE HEADER FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#ifndef FluxLimitedFV_H_
#define FluxLimitedFV_H_

#include <Eigen/Core>
#include <iostream>
#include <stdexcept>

namespace FluxLimitedFV {

/* SAM_LISTING_BEGIN_1 */
template <typename FLUXLIM = std::function<double(double)>>
Eigen::VectorXd fluxlimAdvection(
    double beta, const Eigen::VectorXd& mu0, double h, double tau,
    unsigned int nb_timesteps,
    FLUXLIM&& phi = [](double /*theta*/) { return 1.0; }) {
  if (beta < 0) {
    throw std::domain_error("fluxlimAdvection: negative beta!");
  }
  Eigen::VectorXd mu;  // return vector

  // The length of the vector passing the initial conditions also determines teh
  // number of spatial (dual) cells
  int N = mu0.size();
  // Ratio of spatial and temporal timestep
  double gamma = tau / h;
  // Set initial conditions
  mu = mu0;
#if SOLUTION
  // Constant continuation index mapping tool. This is only one of many ways how
  // to implement the "endpoint conditions"
  auto idx_map = [N](int idx) {
    if (idx < 0) {
      return 0;
    }
    if (idx > N - 1) {
      return N - 1;
    }
    return idx;
  };

  Eigen::VectorXd mu_next(N);
  for (int k = 0; k < nb_timesteps; k++) {
    // Computing forward step
    for (int j = 0; j < N; j++) {
      // Implementation of \prbeqref{eq:ffadv}.
      const double backward_diff = mu(idx_map(j)) - mu(idx_map(j - 1));
      const double shifted_backward_diff =
          mu(idx_map(j - 1)) - mu(idx_map(j - 2));
      const double forward_diff = mu(idx_map(j + 1)) - mu(idx_map(j));
      // Values of flux limiter
      double phi_1, phi_2;
      if (forward_diff == 0.0) {
        phi_1 = 0.0;
      } else {
        phi_1 = phi(backward_diff / forward_diff);
      }
      if (backward_diff == 0.0) {
        phi_2 = 0.0;
      } else {
        phi_2 = phi(shifted_backward_diff / backward_diff);
      }
      // Main recursion
      mu_next(j) = mu(idx_map(j)) - beta * gamma * backward_diff -
                   0.5 * beta * gamma * (1 - beta * gamma) *
                       (phi_1 * forward_diff - phi_2 * backward_diff);
    }
    mu = mu_next;  // step forward (loop update)
  }
#else
  // ========================================
  // Solution code goes here
  // ========================================
#endif
  return mu;
};
/* SAM_LISTING_END_1 */

// Evaluation of quotient of state differences avoiding division by zero
inline double thetaquotient(double u, double v, double w) {
  const double denom = w - v;
  return (denom == 0.0) ? 1E17 * (v - u) : (v - u) / denom;
}

/* SAM_LISTING_BEGIN_2 */
template <typename FLUXLIM = std::function<double(double)>>
Eigen::VectorXd fluxlimBurgers(
    const Eigen::VectorXd& mu0, double h, double tau, unsigned int nb_timesteps,
    FLUXLIM&& phi = [](double /*theta*/) { return 1.0; }) {
  Eigen::VectorXd mu;  // return vector
  int N = mu0.size();  // Number of sptial dual cells
  double gamma = tau / h;
  // Set initial conditions
  mu = mu0;
  // Flux function for Burgers equation
  auto f = [](double u) { return 0.5 * u * u; };
  // Godunov numerical flux
  auto godnfnburgers = [f](double v, double w) -> double {
#if SOLUTION
    if (v > w) {
      if (v + w > 0) {
        return f(v);
      }
      return f(w);
    }
    if (v > 0) {
      return f(v);
    }
    if (0 < w) {
      return 0.0;
    }
    return f(w);

#else
    // ========================================
    // Solution code goes here
    // Replace the following dummy return value:
    return 0.0;
    // ========================================
#endif
  };

#if SOLUTION
  // Constant continuation index mapping tool
  auto idx_map = [N](int idx) {
    if (idx < 0) {
      return 0;
    }
    if (idx > N - 1) {
      return N - 1;
    }
    return idx;
  };
  // Rankine-Hugoniot speed
  auto s_dot = [f](double v, double w) {
    if (v == w) {
      return 0.0;
    }
    return (f(w) - f(v)) / (w - v);
  };
  // The quantity $\theta$ from \prbeqref{eq:theta}
  // Implementation avoid division by zero, see thetaquotient()
  auto theta = [s_dot](double x, double v, double w, double y) {
    if (s_dot(v, w) < 0) {
      if (y == w) {
        return 1.0e17 * (w - v);
      }
      return (w - v) / (y - w);
    }
    if (v == w) {
      return 1.0e17 * (v - x);
    }
    return (v - x) / (w - v);
    ;
  };
  // See \prbeqref{eq:fho}
  auto f_D = [&](int i) -> double {
    double diff = mu(idx_map(i + 1)) - mu(idx_map(i));
    double s_dot_abs = std::abs(s_dot(mu(idx_map(i)), mu(idx_map(i + 1))));
    return 0.5 * s_dot_abs * (1 - gamma * s_dot_abs) * diff *
           phi(theta(mu(idx_map(i - 1)), mu(idx_map(i)), mu(idx_map(i + 1)),
                     mu(idx_map(i + 2))));
  };
  // Auxility arrays: must be declared outside the main timestepping loop
  Eigen::VectorXd mu_next(N);
  Eigen::VectorXd fluxdiff(N);
  for (int k = 0; k < nb_timesteps; k++) {
    // Computing modified flux difference vector, see \prbeqref{eq:flfvf}
    for (int i = 0; i < N; i++) {
      fluxdiff(i) = godnfnburgers(mu(idx_map(i)), mu(idx_map(i + 1))) + f_D(i) -
                    godnfnburgers(mu(idx_map(i - 1)), mu(idx_map(i))) -
                    f_D(i - 1);
    }
    // Computing forward Euler step \prbeqref{eq:ff}
    for (int j = 0; j < N; j++) {
      mu_next(j) = mu(j) - gamma * fluxdiff(j);
    }
    mu = mu_next;  // step forward (loop update)
  }
#else
  // ========================================
  // Solution code goes here
  // ========================================
#endif
  return mu;
}
/* SAM_LISTING_END_2 */

}  // namespace FluxLimitedFV

#endif
