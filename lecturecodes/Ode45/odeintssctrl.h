///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_0 */
// Auxiliary function: default norm for an \eigen vector type
template <class State>
double _norm(const State& y) {
  return y.norm();
}
// Adaptive single-step integrator
template <class DiscEvolOp, class State,
          class NormFunc = decltype(_norm<State>)>
std::vector<std::pair<double, State>> odeintssctrl(
    DiscEvolOp&& Psilow, unsigned int p, DiscEvolOp&& Psihigh, const State& y0,
    double T, double h0, double reltol, double abstol, double hmin,
    NormFunc& norm = _norm<State>) {
  double t = 0;   // initial time $\cob{t_0=0}$\Label[line]{odeintadapt:1}
  State y = y0;   // current state, initialized here
  double h = h0;  // timestep to start with
  std::vector<std::pair<double, State>>
      states;  // vector $\cob{\left(t_k,\Vy_k\right)_k}$
  states.push_back({t, y});

  // Main timestepping loop
  while ((states.back().first < T) && (h >= hmin)) {  // \Label[line]{ssctrl:2}
    State yh =
        Psihigh(h, y);  // high order discrete evolution
                        // \Blue{$\widetilde{\Psibf}^h$}\Label[line]{ssctrl:3}
    State yH = Psilow(h, y);  // low order discrete evolution
                              // \Blue{${\Psibf}^h$}\Label[line]{ssctrl:4}
    double est = norm(
        yH -
        yh);  // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\Label[line]{ssctrl:5}
    double tol =
        std::max(reltol * norm(y),
                 abstol);  // effective tolerance \Label[line]{ssctrl:6a}

    // Optimal stepsize according to \eqref{eq:ssc}
    if (est < tol) {  // step \Magenta{accepted}
                      // \Label[line]{ssctrl:7}\Label[line]{ssctrl:6}
      states.push_back({t = t + std::min(T - t, h),
                        y = yh});  // store next approximate state
    }
    h *= std::max(
        0.5,
        std::min(2., 0.9 * std::pow(tol / est,
                                    1. / (p + 1))));  // \Label[line]{ssctrl:6b}
    if (h < hmin) {
      std::cerr
          << "Warning: Failure at t=" << states.back().first
          << ". Unable to meet integration tolerances without reducing the step"
          << " size below the smallest value allowed (" << hmin
          << ") at time t." << std::endl;
    }
  }
  return states;
}
/* SAM_LISTING_END_0 */
