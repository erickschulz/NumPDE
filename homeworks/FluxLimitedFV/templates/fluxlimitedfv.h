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
  // ========================================
  // Solution code goes here
  // ========================================
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
    // ========================================
    // Solution code goes here
    // Replace the following dummy return value:
    return 0.0;
    // ========================================
  };

  // ========================================
  // Solution code goes here
  // ========================================
  return mu;
}
/* SAM_LISTING_END_2 */

}  // namespace FluxLimitedFV

#endif
