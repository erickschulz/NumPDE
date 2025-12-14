/// Copyright (c) 2016 NumCSE @ ETH Zürich
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

//! \file ode45.hpp Contains header only class for adaptive 4-5 Runge Kutta
//! integration.

/*** BEGIN AUXILIARY FUNCTIONS ***/

//! \brief Compute the norm for a Eigen compatible vector.
//! Wrapper around Eigen norm method. Returns lpNorm<Eigen::Infinity>
//! (the \f$ l^\infty \f$ norm of an Eigen::VectorXx).
//! \tparam T Must have a member method lpNorm<Eigen::Infinity>.
//! \param[in] t The vector you wish to compute the norm of.
//! \return The \f$ l^\infty \f$-norm of t.
template <class T,
          typename std::enable_if<!std::numeric_limits<T>::is_specialized,
                                  bool>::type = true>
inline typename T::Scalar norm_t(const T& t) {
  return t.template lpNorm<Eigen::Infinity>();
}

//! \brief Compute the norm for a scalar.
//! Wrapper around std::abs norm method.
//! \tparam T Must have a T std::abs(T) overload.
//! \param[in] t The scalar you wish to compute the norm of.
//! \return The \f$ l^\infty \f$-norm of t.
template <class T,
          typename std::enable_if<std::numeric_limits<T>::is_specialized,
                                  bool>::type = false>
inline T norm_t(const T& t) {
  return std::abs(t);
}

class termination_error : public std::exception {
  virtual const char* what() const throw() {
    return "Integration terminated prematurely.";
  }
};

/*** END AUXILIARY FUNCTIONS ***/

//! \brief Class for Runge-Kutta-Fehlberg numerical integration.
//! based on nested explicit Runge-Kutta methods of order 4 and 5
//! This class is meant to emulate MATLAB's ode45 integrator.
//! Computes the solution of an ODE with an adaptive Runge Kutta method of
//! order 4/5, for the IVP
//! \f$ y' = f(y), y(0) = y0 \f$, with error control using the difference
//! between the RK4 and RK5 methods.
//!
//! The class works as follows:
//!
//! 1. Call constructor templated by state type and taking the right-hand-side
//! vector field as argument, for example:
//!     Ode45<Eigen::VectorXd> O(f);
//!
//! 2. (optional) Set options, for instance:
//!     O.options.<name> = <value>;
//!
//! 3. Solve phase, e.g.:
//!     auto sol = O.solve(y0, T);
//! or
//!     std::vector<std::pair<RhsType, double>> sol = O.solve(y0, T);
//! in addition a norm can be passed as third argument.
//!
//! 4. (optional) Get statistics:
//!     O.statistics.<stat_you_want_to_get>
//!
//! 5. (optional) Print data (statistics and options):
//!     O.print();
//!
//! \tparam StateType type of the initial data (\f$y_0\f$) and of the
//! solution \f$y(t)\f$. We require to be able to
//! do basic operations on StateSpace, like copy/copy-construct. Furthermore, it
//! is assumed that we have a vector space structure implemented
//! using operators +, *, +=, *=. Moreover we require a norm StateType.norm_t().
//! \tparam RhsType type of the r.h.s. function \f$f\f$, providing
//!      StateType operator()(const StateType & y)
//!
//! The actual embedded Runge-Kutta-Fehlberg method can be selected by the
//! preprocessor flag MATLABCOEFF. If set, uses MATLAB's integrator with 7
//! internal stages. If not set uses a 6-stage embedded method by Petzold &
//! Asher
template <class StateType,
          class RhsType = std::function<StateType(const StateType&)>>
class Ode45 {
 public:
  //! \brief Initialize the class by providing a r.h.s.
  //! RhsType is automatically deduced by the constructor at
  //! construction time.
  //! Copy of r.h.s (of \f$ y'(t) = rhs((y(t)) \f$) is stored internally.
  //! \param[in] f function for the computation of r.h.s.
  //! (e.g. a lambda function).
  Ode45(const RhsType& rhs) : f(rhs), t(0.0) { /* EMPTY */ }

  //! \brief Performs solutions of IVP up to specified final time.
  //! Evolves ODE with initial data \f$y0\f$, up to time \f$T\f$ or
  //! until the ODE integrator breaks down.
  //! \tparam NormFunc Function type for norm function.
  //! \param[in] y0 initial data \f$y_0 = y(0)\f$.
  //! \param[in] T final time for the integration (initial time = 0)
  //! \param[in] norm optional norm function (if a custom norm is needed or
  //! norm_t is not defined, i.e. we use a custom vector type).
  //! \return vector of pairs \f$ (y(t), t) \f$ at snapshot times.
  template <class NormFunc = decltype(norm_t<StateType>)>
  std::vector<std::pair<StateType, double>> solve(
      const StateType& y0, double T, const NormFunc& norm = norm_t<StateType>);

  //! \brief Print statistics and options of this class instance.
  void print();

  //! \brief Stores configuration parameters (a.k.a. options).
  //! Setting values here configures the Ode45 class to use the selected
  //! options.
  //! Setup this before calling solve(...), e.g. with
  //!     Ode45<StateType> O(f);
  //!     Ode45.options.rtol = 10e-5;
  struct Options {
    //!< Set true if you want to save the initial data
    bool save_init = true;
    //!< TODO: Set true if you want a fixed step size
    bool fixed_stepsize = false;
    //!< Set the maximum number of rejected iterations
    unsigned int max_iterations = 5000;
    //!< Set the minimum step size (-1 for none)
    double min_dt = -1.;
    //!< Set the maximum step size (-1 for none)
    double max_dt = -1.;
    //!< Set an initial step size
    double initial_dt = -1.;
    //!< Set a starting time
    double start_time = 0;
    //!< Relative tolerance for the error.
    double rtol = 1e-6;
    //!< Absolute tolerance for the error.
    double atol = 1e-8;
    //!< Set to true before solving to save statistics
    bool do_statistics = false;
    //!< TODO: Perform runtime measurements.
    bool do_timings = false;
    //!< Print more output.
    bool verbose = false;
  } options;

  //! \brief Contain usage statistics.
  //! This will be written (i.e. contain meaningful values)
  //! after a call of solve(...), if do_statistics is set to true.
  struct Statistics {
    //!< Number of loops (sum of all accepted and rejected steps)
    unsigned int cycles = 0;
    //!< Number of actual time steps performed (accepted step)
    unsigned int steps = 0;
    //!< Number of rejected steps per step
    unsigned int rejected_steps = 0;
    //!< Function calls
    unsigned int funcalls = 0;
  } statistics;

 private:
  // A copy of rhs stored during initialization
  RhsType f;
  // Current time
  double t;

  // RK45 coefficients and data, cf. https://github.com/rngantner/
  ////////////////////////////////////////
  //// 20071016, reported by Luis Randez
  // The Runge-Kutta-Fehlberg 4(5) coefficients
  // Coefficients proved on 20060827
  // See p.91 in Ascher & Petzold
  ////////////////////////////////////////
  // Power factor \Blue{$\frac{1}{p+1}$}for error control, \Blue{$p$} = order of
  // lower order methpd
  static constexpr double _pow = 1. / 5;
  // Number of stages
#ifdef MATLABCOEFF
  static const unsigned int _s = 7;
#else
  static const unsigned int _s = 6;
#endif
  // Matrix \Blue{$\FA$} from the Butcher scheme
  static Eigen::MatrixXd _mA;
  // Quadrature weights b_i
  static Eigen::VectorXd _vb4, _vb5;
  // Coefficients c_i relevant for non-autonomous ODEs
  static Eigen::VectorXd _vc;
};

// Matrix \Blue{$\FA$} in Butcher scheme \eqref{eq:BSexpl}
template <class StateType, class RhsType>
Eigen::MatrixXd Ode45<StateType, RhsType>::_mA =
    (Eigen::MatrixXd(Ode45<StateType, RhsType>::_s,
                     Ode45<StateType, RhsType>::_s - 1)
         <<
#ifdef MATLABCOEFF
         0,
     0, 0, 0, 0, 0, 1. / 5., 0, 0, 0, 0, 0, 3. / 40., 9. / 40., 0, 0, 0, 0,
     44. / 45., -56. / 15., 32. / 9., 0, 0, 0, 19372. / 6561., -25360. / 2187.,
     64448. / 6561., -212. / 729., 0, 0, 9017. / 3168., -355 / 33,
     46732. / 5247., 49. / 176., -5103. / 18656., 0, 35. / 384., 0,
     500. / 1113., 125. / 192., -2187. / 6784., 11. / 84.
#else
         0,
     0, 0, 0, 0, 1. / 4., 0, 0, 0, 0, 3. / 32., 9. / 32., 0, 0, 0, 1932. / 2197,
     -7200. / 2197, 7296. / 2197, 0, 0, 439. / 216, -8, 3680. / 513,
     -845. / 4104, 0, -8. / 27., 2, -3544. / 2565, 1859. / 4104, -11. / 40.
#endif
     )
        .finished();

// Quadrature weights \Blue{$b_i$} for the 4th-order method
template <class StateType, class RhsType>
Eigen::VectorXd Ode45<StateType, RhsType>::_vb4 =
    (Eigen::VectorXd(Ode45<StateType, RhsType>::_s) <<
#ifdef MATLABCOEFF
         5179. / 57600.,
     0, 7571. / 16695., 393. / 640., -92097. / 339200., 187. / 2100., 1. / 40.
#else
         25. / 216,
     0, 1408. / 2565, 2197. / 4104, -1. / 5, 0
#endif
     )
        .finished();

// Quadrature weights \Blue{$b_i$} for the 5th-order method
template <class StateType, class RhsType>
Eigen::VectorXd Ode45<StateType, RhsType>::_vb5 =
    (Eigen::VectorXd(Ode45<StateType, RhsType>::_s) <<
#ifdef MATLABCOEFF
         35. / 384.,
     0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0
#else
         16. / 135,
     0, 6656. / 12825, 28561. / 56430, -9. / 50, 2. / 55
#endif
     )
        .finished();

// The coefficients \Blue{$c_i$}, relevant for non-autonomous ODEs.
// Can be computed via row sums of \Blue{$FA$}
template <class StateType, class RhsType>
Eigen::VectorXd Ode45<StateType, RhsType>::_vc =
    Ode45<StateType, RhsType>::_mA.rowwise().sum();

// Stage number of the RK scheme
template <class StateType, class RhsType>
const unsigned int Ode45<StateType, RhsType>::_s;

// solve(): main timestepping method, for autonomous ODEs only
template <class StateType, class RhsType>
template <class NormFunc>
std::vector<std::pair<StateType, double>> Ode45<StateType, RhsType>::solve(
    const StateType& y0, double T, const NormFunc& norm) {
  const double epsilon = std::numeric_limits<double>::epsilon();
  // Setup step size default values if not provided by user
  t = options.start_time;
  unsigned int default_nsteps = 100;
  unsigned int default_minsteps = 10;
  if (options.initial_dt == -1.) {
    options.initial_dt = (T - t) / default_nsteps;
  }
  if (options.max_dt == -1.) {
    options.max_dt = (T - t) / default_minsteps;
  }
  if (options.min_dt == -1.) {
    options.min_dt = (T - t) * epsilon;
  }

  // Vector for returning solution \Blue{$(t_k,\Vy_k)$}
  std::vector<std::pair<StateType, double>> snapshots;

  // The desired (initial) timestep size
  double dt = options.initial_dt;
  if (dt <= 0) {
    std::stringstream ss;
    ss << "Invalid option, dt must be positive (was " << dt << ")!";
    throw std::invalid_argument(ss.str());
  }

  // Push initial data
  if (options.save_init) {
    snapshots.push_back(std::make_pair(y0, t));
  }

  // Temporary containers
  StateType ytemp0 = y0;
  StateType ytemp1 = y0;
  StateType ytemp2 = y0;
  // Pointers for swapping of temporary containers
  StateType *yprev = &ytemp0, *y4 = &ytemp1, *y5 = &ytemp2;

  // Increments \Blue{$\Vk_i$}
  std::vector<StateType> mK;
  mK.resize(_s);

  // Usage statistics
  unsigned int iterations = 0;  // Iterations for current step

  // Main loop, exit if dt too small or final time reached
  while (t < T && dt >= options.min_dt) {
    // Force hitting the endpoint of the time slot exactly
    if (t + dt > T) dt = T - t;
    // Compute the Runge-Kutta increments using the
    // coefficients provided in _mA, _vb, _vc
    mK.front() = f(*yprev);
    for (unsigned int j = 1; j < _s; ++j) {
      mK.at(j) = *yprev;
      for (unsigned int i = 0; i < j; ++i) {
        mK.at(j) += (dt * _mA(j, i)) * mK.at(i);
      }
      mK.at(j) = f(mK.at(j));
    }
    statistics.funcalls += _s;

    // Compute the 4th and the 5th order approximations
    *y4 = *yprev;
    *y5 = *yprev;
    for (unsigned int i = 0; i < _s; ++i) {
      *y4 += (dt * _vb4(i)) * mK.at(i);
      *y5 += (dt * _vb5(i)) * mK.at(i);
    }

    double tau = 2., delta = 1.;
    // Calculate the absolute local truncation error and the acceptable  error
    if (!options.fixed_stepsize) {  // if (!fixed_stepsize)
      delta =
          norm(*y5 - *y4);  // estimated 1-step error \Blue{$\mathtt{EST}_k$}
      tau = std::max(options.rtol * norm(*yprev), options.atol);
    }

    // Check if step is \com{accepted}, if so, advance
    if (delta <= tau) {
      t += dt;
      snapshots.push_back(std::make_pair(*y5, t));
      std::swap(y5, yprev);
      ++statistics.steps;
      statistics.rejected_steps += iterations;
      iterations = 0;
    }

    // Update the step size for the next integration step
    if (!options.fixed_stepsize) {
      if (delta <= std::numeric_limits<double>::epsilon()) {
        dt *= 2;
      } else {
        // Apply formula \eqref{eq:ssc}
        dt *= 0.8 * std::pow(tau / delta, _pow);
      }
      dt = std::min(options.max_dt, dt);
    }
    ++iterations;
    ++statistics.cycles;

    // Check if maximum number of iterations have been exeded
    if (iterations >= options.max_iterations) {
      std::cerr << "Fatal error: the solver has not been successful. The"
                << " integration loop exited at time t = " << t
                << " before the endpoint at tend = " << T
                << " was reached. This happened because the "
                << " maximum number of iterations " << options.max_iterations
                << " was reached. Try to reduce the value of"
                << " \"initial_dt\" and/or \"max_dt\", or increase the"
                << " value of \"max_iterations\". Your ODE may be ill-posed."
                << std::endl;
      throw termination_error();
    }
  }

  // Check if there was a premature exit
  if (t < T) {
    std::cerr << "Fatal error: the solver exited prematurely."
              << " The integration loop exited at time t = " << t
              << " before the endpoint at tend = " << T
              << " was reached. This may happen if the step size becomes"
              << " smaller than the size defined in \"min_dt\"."
              << " Try to reduce the value of "
              << " \"initial_dt\" and/or \"max_dt\"." << std::endl;
    throw termination_error();
  }

  // Returns all collected snapshots
  // TODO: option to select which snapshot to save
  return snapshots;
}

template <class StateType, class RhsType>
void Ode45<StateType, RhsType>::print(void) {
  std::cout << "----------------------------------" << std::endl;
  std::cout << "--- Report of ODE solve Ode45. ---" << std::endl;
#ifdef MATLABCOEFF
  std::cout << "--- MATLAB's Runge-Kutta-Fehlberg method used ---" << std::endl;
#else
  std::cout << "--- BOOST Runge-Kutta-Fehlberg method used ----" << std::endl;
#endif
  std::cout << "----------------------------------" << std::endl;
  std::cout << " + Data:" << std::endl;
  std::cout << "    - current time of simulation:         " << t << std::endl;
  std::cout << " + Options:" << std::endl;
  std::cout << "    - relative tolerance:                 " << options.rtol
            << std::endl;
  std::cout << "    - absolute tolerance:                 " << options.atol
            << std::endl;
  std::cout << "    - minimal stepsize:                   " << options.min_dt
            << std::endl;
  std::cout << "    - maximal stepsize:                   " << options.max_dt
            << std::endl;
  std::cout << "    - initial stepsize:                   "
            << options.initial_dt << std::endl;
  std::cout << "    - max. allowed rejected steps:        "
            << options.max_iterations << std::endl;
  std::cout << "    - starting time:                      "
            << options.start_time << std::endl;
  std::cout << "    - save initial data y(0):             " << options.save_init
            << std::endl;
  std::cout << "    - use fixed stepsize:                 "
            << options.fixed_stepsize << std::endl;
  if (!options.do_statistics) return;
  std::cout << " + Statistics:" << std::endl;
  std::cout << "    - number of steps:                    " << statistics.steps
            << std::endl;
  std::cout << "    - number of rejected steps:           "
            << statistics.rejected_steps << std::endl;
  std::cout << "    - number of (while-)loops:            " << statistics.cycles
            << std::endl;
  std::cout << "    - function calls:                     "
            << statistics.funcalls << std::endl;
  std::cout << " + Timing:" << std::endl;
  std::cout << "----------------------------------" << std::endl;
  std::cout << "----------------------------------" << std::endl;
  // TODO: time solver
}
