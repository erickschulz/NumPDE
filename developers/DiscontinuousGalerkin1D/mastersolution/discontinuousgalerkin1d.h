/**
 * @file discontinuousgalerkin1d.h
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <iostream>
#include <utility>

namespace DiscontinuousGalerkin1D {

/**
 * @brief Returns the mass matrix B according to the exercise.
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return diagonal square matrix of dimension 2 * (Ml + Mr + 1)
 */
Eigen::SparseMatrix<double> compBmat(int Ml, int Mr, double h);

/**
 * @brief Computes G according to the problem description.
 * @param mu expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * @param f flux function matching std::function<double(double)>
 * @param F numerical flux matching std::function<double(double, double)>
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return vector of length 2 * (Ml + Mr + 1) approximating G(mu)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR, typename NUMFLUX>
Eigen::VectorXd G(const Eigen::VectorXd& mu, FUNCTOR&& f, NUMFLUX&& F, int Ml,
                  int Mr, double h) {
  const int N_half = (Ml + Mr + 1);
  const int N = 2 * N_half;
  Eigen::VectorXd Gvec(N);
#if SOLUTION
  double uN_xminus = 0.0;  // since we extend mu to the left by zero
  double uN_xplus = mu(0) - 0.5 * h * mu(1);
  double F_old;
  double F_new = F(uN_xminus, uN_xplus);
  const double w = h / (2.0 * std::sqrt(3.0));

  for (int i = 0; i < N_half - 1; ++i) {
    uN_xminus = mu(2 * i) + 0.5 * h * mu(2 * i + 1);
    uN_xplus = mu(2 * (i + 1)) - 0.5 * h * mu(2 * (i + 1) + 1);
    F_old = F_new;
    F_new = F(uN_xminus, uN_xplus);
    Gvec(2 * i) = F_new - F_old;

    double x_minus = mu(2 * i) - w * mu(2 * i + 1);
    double x_plus = mu(2 * i) + w * mu(2 * i + 1);
    double I = 0.5 * h * (f(x_minus) + f(x_plus));
    Gvec(2 * i + 1) = 0.5 * h * (F_new + F_old) - I;
  }

  uN_xminus = mu(2 * (N_half - 1)) + 0.5 * h * mu(2 * (N_half - 1) + 1);
  uN_xplus = 0.0;  // since we extend mu to the right by zero
  F_old = F_new;
  F_new = F(uN_xminus, uN_xplus);
  Gvec(2 * (N_half - 1)) = F_new - F_old;

  double x_minus = mu(2 * (N_half - 1)) - w * mu(2 * (N_half - 1) + 1);
  double x_plus = mu(2 * (N_half - 1)) + w * mu(2 * (N_half - 1) + 1);
  double I = 0.5 * h * (f(x_minus) + f(x_plus));
  Gvec(2 * (N_half - 1) + 1) = 0.5 * h * (F_new + F_old) - I;
#else
  //====================
  // Your code goes here
  // Fill the vector Gvec
  //====================
#endif
  return Gvec;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Time evolution by an explicit 2nd order Runge-Kutta method.
 * @param mu0 expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * representing the initial data
 * @param f flux function matching std::function<double(double)>
 * @param F numerical flux matching std::function<double(double, double)>
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return expansion coefficients vector of length 2 * (Ml + Mr + 1)
 * representing the solution at time T
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR, typename NUMFLUX>
Eigen::VectorXd dgcl(Eigen::VectorXd mu0, FUNCTOR&& f, NUMFLUX&& F, double T,
                     int Ml, int Mr, double h, unsigned int m) {
#if SOLUTION
  Eigen::SparseMatrix<double> B = compBmat(Ml, Mr, h);
  Eigen::SparseMatrix<double> Binv = B.cwiseInverse();

  auto G_bound = [&f, &F, T, Ml, Mr, h](const Eigen::VectorXd& mu) {
    return G(mu, std::forward<FUNCTOR>(f), std::forward<NUMFLUX>(F), Ml, Mr, h);
  };
  // Timestepping based on explicit midpoint rule, a 2-stage explicit Rune-Kutta
  // single-step method
  double tau = T / m;
  for (int i = 0; i < m; ++i) {
    // First compute the increment and then update the state vector, see
    // \lref{eq:Chemp}
    Eigen::VectorXd k = -Binv * G_bound(mu0);
    mu0 = mu0 - tau * Binv * G_bound(mu0 + 0.5 * tau * k);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return mu0;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Engquist-Osher numerical flux for the flux function f(u) = u * (1 -
 * u).
 * @param mu0 expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * representing the initial data
 * @param v left state
 * @param w right state
 * @return F(v, w)
 */
double Feo(double v, double w);

struct Solution {
  Solution(const Solution& other) : x_(other.x_), u_(other.u_) {
    std::cout << "Called copy contructor" << std::endl;
  }
  Solution(Eigen::VectorXd x, Eigen::VectorXd u)
      : x_(std::move(x)), u_(std::move(u)) {}
  Eigen::VectorXd x_;
  Eigen::VectorXd u_;
};

/**
 * @brief Time evolution according to dgcl(...) on the spacial interval [-2, 2]
 * with mesh-width h = 0.05, on th time interval [0, 1] with timestep size h
 * / 3. The solution at endtime is written to solution.csv.
 */
Solution solveTrafficFlow();

}  // namespace DiscontinuousGalerkin1D
