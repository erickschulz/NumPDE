///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_0 */
// Auxiliary function: default norm for an \eigen vector type
template <class State>
double _norm(const State& y) {
  return y.norm();
}

// Adaptive numerical integrator based on local-in-time stepsize control
template <class DiscEvolOp, class State,
          class NormFunc = decltype(_norm<State>)>
std::vector<std::pair<double, State>> odeintadapt(
    DiscEvolOp&& Psilow, DiscEvolOp&& Psihigh, const State& y0, double T,
    double h0, double reltol, double abstol, double hmin,
    NormFunc& norm = _norm<State>) {
  double t = 0;   // initial time $\cob{t_0=0}$\Label[line]{odeintadapt:1}
  State y = y0;   // current state
  double h = h0;  // timestep to start with
  std::vector<std::pair<double, State>>
      states;                // vector of times/computed states:
                             // $\cob{\left(t_k,\Vy_k\right)_k}$
  states.push_back({t, y});  // initial time and state

  while ((states.back().first < T) &&
         (h >= hmin)) {  // \Label[line]{odeintadapt:2}
    State yh = Psihigh(
        h, y);  // high order discrete evolution \Blue{$\widetilde{\Psibf}^h$}
                // \Label[line]{odeintadapt:3}
    State yH = Psilow(h, y);  // low order discrete evolution
                              // \Blue{${\Psibf}^h$} \Label[line]{odeintadapt:4}
    double est =
        norm(yH - yh);  // local error estimate
                        // \Blue{$\mathrm{EST}_k$}\Label[line]{odeintadapt:5}

    if (est <
        std::max(
            reltol * norm(y),
            abstol)) {  // step \Magenta{accepted} \Label[line]{odeintadapt:6}
      y = yh;           // use high order approximation
      t = t + std::min(T - t, h);  // next time \Blue{$t_k$}
      states.push_back({t, y});    // \Label[line]{odeintadapt:7}
      h = 1.1 * h;  // try with increased stepsize \Label[line]{odeintadapt:8}
    } else {        // step \Magenta{rejected}
      h = h / 2;    // try with half the stepsize \Label[line]{odeintadapt:9}
    }
    // Numerical integration has ground to a halt !
    if (h < hmin) {
      std::cerr << "Warning: Failure at t=" << states.back().first
                << ". Unable to meet integration tolerances without reducing "
                   "the step "
                << "size below the smallest value allowed (" << hmin
                << ") at time t." << std::endl;
    }
  }
  return states;
}
/* SAM_LISTING_END_0 */
