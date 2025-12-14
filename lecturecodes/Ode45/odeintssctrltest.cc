#include <cmath>
#include <functional>
#include <iostream>

#include "odeintssctrl.h"

int main() {
  // A scalar initial-value problem
  using State_t = double;
  using DiscEvolOp = std::function<State_t(double, State_t)>;
  // Definition of right-hand-side for ODE
  auto f = [](double y) { return y * y; };
  // Initial value
  double y0 = 0.5;
  // Exact solution
  auto y = [](double t) { return 1.0 / (2.0 - t); };
  // Final time
  double T = 1.9;
  // The norm is just the modulus in the scalar case
  auto norm = [](double x) { return std::fabs(x); };

  // Low-order method: explicit euler (order 1)
  DiscEvolOp psilow = [&](double h, double y) { return y + h * f(y); };

  // "High-order method": explicit trapezoidal rule (order 2)
  DiscEvolOp psihigh = [&](double h, double y) {
    double k1 = f(y);
    double k2 = f(y + h * k1);
    return y + (h / 2.) * (k1 + k2);
  };
  // Call simple adaptive integrator
  std::vector<std::pair<double, double>> states =
      odeintssctrl(psilow, 1, psihigh, y0, 1.9, 0.2, 1e-3, 1e-4, 1e-4, norm);
  // Output solution and error
  std::cout << "Adaptive integration, " << states.size() - 1 << " timesteps"
            << std::endl;
  for (auto& ty : states) {
    std::cout << "t = " << ty.first << ": y = " << ty.second
              << ", error = " << norm(ty.second - y(ty.first)) << std::endl;
  }
  return 0;
}
