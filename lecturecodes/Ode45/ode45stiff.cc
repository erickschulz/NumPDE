///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Example {ex:ode45stiff}
/// (C) 2016 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#define MATLABCOEFF true

#include "ode45.h"

void case1() {
  std::cout << "Explicit adaptive RK-SSM for stiff scalar logistic IVP"
            << std::endl;
  /* SAM_LISTING_BEGIN_0 */
  // Types to be used for a scalar ODE with state space \Blue{$\bbR$}
  using StateType = double;
  using RhsType = std::function<StateType(StateType)>;
  // Logistic differential equation \eqref{eq:logode}
  double lambda = 500.0;
  RhsType f = [lambda](StateType y) { return lambda * y * y * (1 - y); };
  StateType y0 = 0.01;  // Initial value, will create a STIFF IVP
  // State space \Blue{$\bbR$}, simple modulus supplies norm
  auto normFunc = [](StateType x) { return fabs(x); };

  // Invoke explicit Runge-Kutta method with stepsize control
  Ode45<StateType, RhsType> integrator(f);
  // Set rather loose tolerances
  integrator.options.rtol = 0.1;
  integrator.options.atol = 0.001;
  integrator.options.min_dt = 1E-18;
  std::vector<std::pair<StateType, double>> states =
      integrator.solve(y0, 1.0, normFunc);
  // Output information accumulation during numerical integration
  integrator.options.do_statistics = true;
  integrator.print();
  /* SAM_LISTING_END_0 */

  for (auto state : states) {
    std::cout << "t = " << state.second << ", y = " << state.first << std::endl;
  }
}

void case2() {
  /* SAM_LISTING_BEGIN_1 */
  // Types to be used for a scalar ODE with state space \Blue{$\bbR$}
  using StateType = double;
  using RhsType = std::function<StateType(StateType)>;
  // Decay differential equation \eqref{eq:logode}
  double lambda = 80.0;
  RhsType f = [lambda](StateType y) { return -lambda * y; };
  StateType y0 = 1.0;  // Initial value, will create a STIFF IVP
  // State space \Blue{$\bbR$}, simple modulus supplies norm
  auto normFunc = [](StateType x) { return fabs(x); };

  // Invoke explicit Runge-Kutta method with stepsize control
  Ode45<StateType, RhsType> integrator(f);
  // Set rather loose tolerances
  integrator.options.rtol = 0.1;
  integrator.options.atol = 0.001;
  integrator.options.min_dt = 1E-18;
  std::vector<std::pair<StateType, double>> states =
      integrator.solve(y0, 1.0, normFunc);
  // Output information accumulation during numerical integration
  integrator.options.do_statistics = true;
  integrator.print();
  /* SAM_LISTING_END_1 */

  for (auto state : states) {
    std::cout << "t = " << state.second << ", y = " << state.first << std::endl;
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << "[1-2]" << std::endl;
    std::cout
        << "Running both cases by default"
        << std::endl;  // This is unconventional but useful for the CI pipeline
    case1();
    case2();
  } else {
    int select;
    try {
      select = std::stoi(argv[1]);
    } catch (const std::exception& e) {
      std::cout << e.what() << std::endl;
      select = 0;
    }
    switch (select) {
      case 1: {
        std::cout << "Explicit adaptive RK-SSM for stiff scalar logistic IVP"
                  << std::endl;
        case1();
        break;
      }
      case 2: {
        std::cout << "Explicit adaptive RK-SSM for stiff decay IVP"
                  << std::endl;
        case2();
        break;
      }
      default: {
        std::cerr << "Usage: " << argv[0] << "[1-2]" << std::endl;
        exit(1);
      }
    }
  }
}
