/**
 * @file laxwendroffscheme.cc
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "laxwendroffscheme.h"

#include <Eigen/Core>
#include <cmath>

namespace LaxWendroffScheme {

namespace Constant {
constexpr double e = 2.71828182845904523536;
constexpr double pi = 3.14159265358979323846;
}  // namespace Constant

constexpr double Square(double x) { return x * x; }

/**
 * @brief Computes the right-hand side according to the Lax-Wendroff scheme.
 * @param mu mu^(k-1) (i.e. mu at timestep k-1)
 * @param gamma tau / h, where tau = timestep size and h = spatial meshwidth
 * @return mu^(k) (i.e. mu at timestep k)
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd LaxWendroffRhs(const Eigen::VectorXd& mu, double gamma) {
  int N = mu.size();
  Eigen::VectorXd result(N);

  auto f = [](double x) { return std::exp(x); };
  auto df2 = [](double x) { return Square(std::exp(x)); };

  // Lax-Wendroff fully discrete evolution \prbeqref{eq:12} for
  // $\cob{f(u)=e^u}$. Store the values of $f(x)$ and $(f'(x))^2$ from the
  // previous iteration:

  // We extend $\texttt{mu}$ to the left by $\texttt{mu(0)}$:
  double mu_left = mu(0);

  // $f'\left(\tfrac{1}{2}\left(\mu_j+\mu_{j-1}\right)\right)^2$:
  double df2_old = df2(0.5 * (mu(0) + mu_left));

  // $f'\left(\tfrac{1}{2}\left(\mu_{j+1}+\mu_j\right)\right)^2$:
  double df2_new = df2(0.5 * (mu(1) + mu(0)));

  double f_old = f(mu_left);  // $f\left(\mu_{j-1}\right)$
  double f_mid = f(mu(0));    // $f\left(\mu_j\right)$
  double f_new = f(mu(1));    // $f\left(\mu_{j+1}\right)$

  result(0) = mu(0) - 0.5 * gamma * (f_new - f_old) +
              0.5 * Square(gamma) *
                  (df2_new * (mu(1) - mu(0)) - df2_old * (mu(0) - mu_left));

  for (int j = 1; j < N - 1; ++j) {
    df2_old = df2_new;
    df2_new = df2(0.5 * (mu(j + 1) + mu(j)));
    f_old = f_mid;
    f_mid = f_new;
    f_new = f(mu(j + 1));
    result(j) =
        mu(j) - 0.5 * gamma * (f_new - f_old) +
        0.5 * Square(gamma) *
            (df2_new * (mu(j + 1) - mu(j)) - df2_old * (mu(j) - mu(j - 1)));
  }

  // We extend $\texttt{mu}$ to the right by $\texttt{mu(N-1)}$:
  double mu_right = mu(N - 1);

  df2_old = df2_new;
  df2_new = df2(0.5 * (mu_right + mu(N - 1)));

  f_old = f_mid;
  f_new = f(mu_right);

  result(N - 1) = mu(N - 1) - 0.5 * gamma * (f_new - f_old) +
                  0.5 * Square(gamma) *
                      (df2_new * (mu_right - mu(N - 1)) -
                       df2_old * (mu(N - 1) - mu(N - 2)));

  return result;
}

Eigen::VectorXd solveLaxWendroff(const Eigen::VectorXd& u0, double T,
                                 unsigned int M) {
  double gamma = 1.0 / Constant::e;
  Eigen::VectorXd mu = u0;
  // Main timestepping loop
  for (int j = 0; j < M; ++j) mu = LaxWendroffRhs(mu, gamma);
  return mu;
}

/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
// Build spatial grid
Eigen::VectorXd getXValues(double T, unsigned int M) {
  double tau = T / M;
  double h = Constant::e * tau;
  int j_max = (int)(std::ceil((3.0 * T + 1.0) / h) + 0.5);
  int j_min = (int)(std::floor(-3.0 * T / h) - 0.5);
  unsigned int N = j_max - j_min + 1;
  return Eigen::VectorXd::LinSpaced(N, j_min * h, j_max * h);
}
Eigen::VectorXd numexpLaxWendroffRP(const Eigen::VectorXi& M) {
  const double T = 1.0;
  const int M_size = M.size();
  Eigen::VectorXd error(M_size);
  // Initial values for the Riemann problem
  auto u_initial = [](double x) { return 0.0 <= x ? 1.0 : 0.0; };
  // Exact solution \prbeqref{eq:solrp} at time $T = 1.0$
  auto u_exact = [](double x) {
    return (x <= 1.0) ? 0.0 : ((Constant::e <= x) ? 1.0 : std::log(x));
  };
  for (int i = 0; i < M_size; ++i) {
    Eigen::VectorXd x = getXValues(T, M(i));
    Eigen::VectorXd u0 = x.unaryExpr(u_initial);
    Eigen::VectorXd uT = solveLaxWendroff(u0, T, M(i));

    double tau = T / M(i);
    double h = Constant::e * tau;
    error(i) = h * (x.unaryExpr(u_exact) - uT).lpNorm<1>();
  }
  return error;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Evaluates the discrete function u at position x by linear
 * interpolation
 * @param u descrete function values at spatial positions y
 * @param y vector of same length as u, representing the nodes of u
 * @return best linear interpolation of u at spacial position x
 */
double eval(const Eigen::VectorXd& u, const Eigen::VectorXd& y, double x) {
  int N = y.size();
  double a = y(0);
  double b = y(N - 1);

  if (x <= a) return u(0);
  if (b <= x) return u(N - 1);

  double lambda = (x - a) / (b - a);
  int k0 = (int)(lambda * (N - 1));
  int k1 = k0 + 1;

  lambda = (x - y(k0)) / (y(k1) - y(k0));
  return lambda * u(k1) + (1.0 - lambda) * u(k0);
}

/* SAM_LISTING_BEGIN_9 */
double smoothU0(double x) {
  return (x < 0.0)
             ? 0.0
             : ((1.0 < x) ? 1.0 : Square(std::sin(0.5 * Constant::pi * x)));
}
Eigen::VectorXd referenceSolution(const Eigen::VectorXd& x) {
  double T = 1.0;
  // Reference solution on a very fine mesh
  unsigned int M = 3200;

  Eigen::VectorXd y = getXValues(T, M);
  Eigen::VectorXd u0 = y.unaryExpr(&smoothU0);
  Eigen::VectorXd u = solveLaxWendroff(u0, T, M);
  int N = x.size();
  Eigen::VectorXd u_ref(N);
  // The vector u is larger than u_ref. Use eval() from above the "evaluate" u
  // at the positions x(i) and thus obtain the reference solution u_ref.
  for (int i = 0; i < N; ++i) {
    u_ref(i) = eval(u, y, x(i));
  }
  return u_ref;
}
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd numexpLaxWendroffSmoothU0(const Eigen::VectorXi& M) {
  const double T = 1.0;
  const int M_size = M.size();
  Eigen::VectorXd error(M_size);

  for (int i = 0; i < M_size; ++i) {
    Eigen::VectorXd x = getXValues(T, M(i));
    Eigen::VectorXd u0 = x.unaryExpr(&smoothU0);
    Eigen::VectorXd uT = solveLaxWendroff(u0, T, M(i));

    double tau = T / M(i);
    double h = Constant::e * tau;
    error(i) = h * (referenceSolution(x) - uT).lpNorm<1>();
  }
  return error;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_7 */
Eigen::VectorXd solveGodunov(const Eigen::VectorXd& u0, double T,
                             unsigned int M) {
  double tau = T / M;
  double h = Constant::e * tau;
  unsigned int N = u0.size();
  Eigen::VectorXd mu = u0;

  for (int i = 0; i < M; ++i) {
    for (int j = N - 1; 0 < j; --j) {
      mu(j) = mu(j) - tau / h * (std::exp(mu(j)) - std::exp(mu(j - 1)));
    }
  }
  return mu;
}

/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
Eigen::VectorXd numexpGodunovSmoothU0(const Eigen::VectorXi& M) {
  const double T = 1.0;
  const int M_size = M.size();
  Eigen::VectorXd error(M_size);

  for (int i = 0; i < M_size; ++i) {
    Eigen::VectorXd x = getXValues(T, M(i));
    Eigen::VectorXd u0 = x.unaryExpr(&smoothU0);
    Eigen::VectorXd uT = solveGodunov(u0, T, M(i));

    double tau = T / M(i);
    double h = Constant::e * tau;
    error(i) = h * (referenceSolution(x) - uT).lpNorm<1>();
  }
  return error;
}
/* SAM_LISTING_END_8 */

}  // namespace LaxWendroffScheme
