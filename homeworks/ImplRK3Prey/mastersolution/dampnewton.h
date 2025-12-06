#ifndef DAMPNEWTON_H_
#define DAMPNEWTON_H_

/**
 * @file dampnewton.h
 * @brief NPDE homework ImplRK3Prey code
 * @author Unknown, Oliver Rietmann
 * @date 29.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/LU>

namespace ImplRK3Prey {

constexpr double min(double x, double y) { return x < y ? x : y; };

template <typename FuncType, typename JacType>
void dampnewton(FuncType&& F, JacType&& DF, Eigen::VectorXd& x,
                double rtol = 1e-4, double atol = 1e-6) {
  const int n = x.size();
  const double lmin = 1E-3;     // Minimal damping factor
  double lambda = 1.0;          // Initial and actual damping factor
  Eigen::VectorXd s(n), st(n);  // Newton corrections
  Eigen::VectorXd xn(n);        // Tentative new iterate
  double sn, stn;               // Norms of Newton corrections

  do {
    auto jacfac = DF(x).lu();  // LU-factorize Jacobian

    s = jacfac.solve(F(x));  // Newton correction
    sn = s.norm();           // Norm of Newton correction
    lambda *= 2.0;
    do {
      lambda /= 2;
      if (lambda < lmin) throw "No convergence: lambda -> 0";
      xn = x - lambda * s;       // {\bf Tentative next iterate}
      st = jacfac.solve(F(xn));  // Simplified Newton correction
      stn = st.norm();
    } while (stn > (1 - lambda / 2) * sn);  // {\bf Natural monotonicity test}
    x = xn;                           // Now: xn accepted as new iterate
    lambda = min(2.0 * lambda, 1.0);  // Try to mitigate damping
  }
  // Termination based on simplified Newton correction
  while ((stn > rtol * x.norm()) && (stn > atol));
}

}  // namespace ImplRK3Prey

#endif  // #define DAMPNEWTON_H_
