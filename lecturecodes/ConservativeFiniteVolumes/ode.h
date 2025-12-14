// Demonstration code for course Numerical Methods for Partial Differential
// Equations Author: R. Hiptmair, O. Rietmann, SAM, ETH Zurich Date: May 2020
#ifndef ODE_HPP
#define ODE_HPP

#include <Eigen/Core>
#include <tuple>
#include <utility>
#include <vector>

// HACK BEGIN
// Avoid compatibility issue between
// boost/numeric/odeint/external/eigen/eigen_algebra.hpp by shaddwing it with
// modified header ./eigen_algebra.hpp
#include "eigen_algebra.hpp"
// HACK END

// Use BOOST numerical library, see
// https://www.boost.org/doc/libs/1_66_0/libs/numeric/odeint/doc/html/index.html
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <tuple>

namespace ConsFV {

/**
 * @brief Tries to mimic ode45 from Matlab
 * Integrates the system of differential equations y' = f(t, y) from time t0 to
 * tfinal with initial conditions y0.
 * @tparam ODEFUN a callable object with signature
 * std::function<void(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const
 * double t)>, which provides the right-hand-side vectorfield f(x,t) for the
 * ODE. The first argument supplies the state x, the second is used for
 * returning the state f(x,t), and the third argument passes the time t.
 * The state is a vector of type Eigen::VectorXd
 *
 * @param odefun The function f(t, y) for a scalar t and a vector y.
 * @param t0 The starting time.
 * @param tfinal The ending time.
 * @param y0 The initial conditions.
 * @param abserr The absolute error tolerance.
 * @param relerr The relative error tolerance.
 * @return temporal grid and sequence of approximate states
 */
template <typename ODEFUN>
std::pair<std::vector<double>, std::vector<Eigen::VectorXd>> ode45(
    ODEFUN&& odefun, const double t0, const double tfinal,
    const Eigen::VectorXd& y0, const double abserr = 1.0E-8,
    const double relerr = 1.0E-6) {
  // initialization
  Eigen::VectorXd mu0 = y0;
  std::vector<double> timesteps{};
  std::vector<Eigen::VectorXd> mus{};
  // initialize a runge kutta stepper with adaptive step size and using the
  // state type Eigen::VectorXd
  auto stepper = boost::numeric::odeint::make_controlled(
      abserr, relerr,
      boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXd>());
  // solve the system recording the time point and corresponding solution vector
  // at every step
  boost::numeric::odeint::integrate_adaptive(
      stepper, odefun, mu0, t0, tfinal, 0.01,
      [&mus, &timesteps](const Eigen::VectorXd& x, double t) {
        mus.resize(mus.size() + 1, x);
        timesteps.resize(timesteps.size() + 1, t);
      });
  return {timesteps, mus};
}

}  // namespace ConsFV

#endif
