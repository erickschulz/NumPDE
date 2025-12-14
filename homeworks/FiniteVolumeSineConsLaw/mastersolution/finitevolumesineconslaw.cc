/**
 * @file finitevolumesineconslaw.cc
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "finitevolumesineconslaw.h"

#include <Eigen/Core>
#include <cmath>

namespace FiniteVolumeSineConsLaw {

/* SAM_LISTING_BEGIN_1 */
constexpr double PI = 3.14159265358979323846;

double f(double x) { return std::sin(PI * x); }

double sineGodFlux(double v, double w) {
  double result;
  // Rankine-Hugoniot speed. No safeguards against cancellation are taken. For v
  // close to w the result will severely be affected by amplified round-off
  // error. However, this does not do any harm, because we are interested in the
  // the sign alone.
  double s = (v != w) ? (f(w) - f(v)) / (w - v) : PI * std::cos(PI * v);
  // Treat all different cases separately
  if (((v < w) && (s > 0.0)) || ((v >= w) && (std::cos(PI * v) > 0.0))) {
    result = f(v);
  } else if (((v < w) && (s < 0.0)) || ((v >= w) && (std::cos(PI * w) < 0.0))) {
    result = f(w);
  } else
    result = f(0.5);

  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd sineClawRhs(const Eigen::VectorXd& mu) {
  int N = mu.size();
  double h = 12.0 / N;
  Eigen::VectorXd result(N);

  double F_minus;
  double F_plus = sineGodFlux(0.0, mu(0));
  for (int j = 0; j < N - 1; ++j) {
    F_minus = F_plus;
    F_plus = sineGodFlux(mu(j), mu(j + 1));
    result(j) = -1.0 / h * (F_plus - F_minus);
  }
  F_minus = F_plus;
  F_plus = sineGodFlux(mu(N - 1), 0.0);
  result(N - 1) = -1.0 / h * (F_plus - F_minus);

  return result;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
bool blowup(const Eigen::VectorXd& mu) {
  return mu.minCoeff() < 0.0 || mu.maxCoeff() > 2.0;
}

unsigned int findTimesteps() {
  const unsigned int N = 600;
  unsigned int ML = 100;
  unsigned int MR = 200;

  for (int i = 0; i < 100; ++i) {
    if (ML == MR) {
      return MR;
    }
    unsigned int M = (unsigned int)(0.5 * (ML + MR));
    if (blowup(solveSineConsLaw(&sineClawRhs, N, M))) {
      ML = M + 1;
    } else {
      MR = M;
    }
  }

  return MR;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd sineClawReactionRhs(const Eigen::VectorXd& mu, double c) {
  Eigen::VectorXd rhs(mu.size());
  rhs = sineClawRhs(mu) - c * mu;
  return rhs;
}
/* SAM_LISTING_END_4 */

}  // namespace FiniteVolumeSineConsLaw
