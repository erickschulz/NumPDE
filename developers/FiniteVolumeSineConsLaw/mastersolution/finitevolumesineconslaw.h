/**
 * @file finitevolumesineconslaw.h
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace FiniteVolumeSineConsLaw {

/**
 * @brief Godunov numerical flux with F, with f(u) = sin(PI * u).
 *
 * @param v 0 <= v <= 1
 * @param w 0 <= w <= 1
 * @return F(v, w)
 */
double sineGodFlux(double v, double w);

/**
 * @brief RHS of the semi-discretized ODE on the spacial interval [-6, 6].
 * The solution u(x, t) is assumed to be zero on the complement.
 *
 * @param mu cell averages
 * @return RHS evaluated at mu
 */
Eigen::VectorXd sineClawRhs(const Eigen::VectorXd& mu);

/**
 * @brief Applies the explicit trapezoidal method to an (semi-discretized) ODE.
 *
 * @param g RHS of the (semi-discretized) ODE
 * @param y0 initial data u(x, 0) for x in [-6, 6]
 * @param T final time T > 0
 * @param M number of time steps
 * @return Cell averages of the approximate solution u(x, T)
 * on the same sub-intervals on which y0 was given.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename RHSFUNCTOR>
Eigen::VectorXd explTrpzTimestepping(RHSFUNCTOR&& g, const Eigen::VectorXd& y0,
                                     double T, unsigned int M) {
  double tau = T / M;
  Eigen::VectorXd y = y0;

#if SOLUTION
  for (int j = 0; j < M; ++j) {
    Eigen::VectorXd k0 = g(y);
    Eigen::VectorXd k1 = g(y + tau * k0);
    y = y + tau * 0.5 * (k0 + k1);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return y;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves the fully discretized ODE by explTrpzTimestepping(...),
 * where the initial data is the characterisitc function of [0, 1),
 * on the spacial interval [-6, 6] and time interval [0, T=1].
 *
 * @param g RHS of the (semi-discretized) ODE
 * @param N number of uniform cells partitioning [-6, 6]
 * @param M number of time steps
 * @return Cell averages of the approximate solution u(x, T=1).
 */
/* SAM_LISTING_BEGIN_2 */
template <typename RHSFUNCTOR>
Eigen::VectorXd solveSineConsLaw(RHSFUNCTOR&& g, unsigned int N,
                                 unsigned int M) {
  const double h = 12.0 / N;
  const double T = 1.0;
  Eigen::VectorXd result(N);

#if SOLUTION
  Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, -6.0 + 0.5 * h, 6.0 - 0.5 * h);
  Eigen::VectorXd u0 =
      x.unaryExpr([](double x) { return 0.0 <= x && x < 1.0 ? 1.0 : 0.0; });
  result = explTrpzTimestepping(std::forward<RHSFUNCTOR>(g), u0, T, M);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return result;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Compute the optimal number of time steps M in the setting above,
 * with N = 600 and T = 1.0.

 * @return Minimal M, sucht that the final solution lies in [0, 2].
 */
unsigned int findTimesteps();

/**
 * @brief Same as sineClawRhs(...), but with the additional reaction term
 * -c * u(x, t).
 *
 * @param mu cell averages
 * @param c c > 0
 * @return RHS evaluated at mu
 */
Eigen::VectorXd sineClawReactionRhs(const Eigen::VectorXd& mu, double c);

}  // namespace FiniteVolumeSineConsLaw
