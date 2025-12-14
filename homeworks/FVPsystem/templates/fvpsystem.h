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

  /* Your code goes here! */

  // Return
  return mu;
}
/* SAM_LISTING_END_1 */

}  // namespace FVPsystem

#endif
