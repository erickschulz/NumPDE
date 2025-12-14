/**
 * @file conslawwithsource.h
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 20.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef CONSLAWWITHSOURCE_H_
#define CONSLAWWITHSOURCE_H_

#include <Eigen/Core>
#include <cmath>
#include <utility>

namespace ConsLawWithSource {

/**
 * @brief Implements the Godunov flux for \Blue{$f(u)=e^u-u$}.
 *
 * @param v left state \Blue{$v$}
 * @param w right state \Blue{$w$}
 * @return value of the Godunov flux \Blue{$F_{GD}(v,w)$}
 */
double godnfn(double v, double w);

/**
 * @brief Implements \Blue{$-h$} times the RHS of the semi-descrete ODE arising
 * in the finite volume discretization in space, in presence of source function.
 *
 * @param mu vector of length N containing the cell averages
 * @param F numerical flux function matching std::function<double(double,
 * double)>
 * @param s source function matching std::function<double(double)>
 * @param h cell size
 * @return vector of length N containing -h times the RHS of the
 * semi-discrete ODE applied to mu
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FunctionF, typename SourceFunction>
Eigen::VectorXd fluxdiffsource(const Eigen::VectorXd& mu, FunctionF&& F,
                               SourceFunction&& s, double h) {
#if SOLUTION
  unsigned n = mu.size();
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);

  // constant continuation to the left by mu[0]
  fd[0] = F(mu[0], mu[1]) - F(mu[0], mu[0]) - h * s(mu[0]);
  for (unsigned j = 1; j < n - 1; ++j) {
    fd[j] = F(mu[j], mu[j + 1]) - F(mu[j - 1], mu[j]) - h * s(mu[j]);
  }
  // constant continuation to the right by mu[n - 1]
  fd[n - 1] =
      F(mu[n - 1], mu[n - 1]) - F(mu[n - 2], mu[n - 1]) - h * s(mu[n - 1]);

  return fd;
#else
  //====================
  // ADAPT THE CODE below so that it handels the source function s correctly
  //====================

  unsigned n = mu.size();                         // length of state vector
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);  // return vector

  // constant continuation of data for \Blue{$x\leq a$}!
  fd[0] = F(mu[0], mu[1]) - F(mu[0], mu[0]);
  for (unsigned j = 1; j < n - 1; ++j) {
    fd[j] = F(mu[j], mu[j + 1]) - F(mu[j - 1], mu[j]);  // see \eqref{eq:2pcf}
  }
  // constant continuation of data for \Blue{$x\geq b$}!
  fd[n - 1] = F(mu[n - 1], mu[n - 1]) - F(mu[n - 2], mu[n - 1]);
  // Efficient thanks to return value optimization (RVO)
  return fd;
#endif
}
/* SAM_LISTING_END_1 */

/**
 * @brief Computes the time evolution of the total mass associated to the
 * solution of the finite volume method, on the time interval [0,3] and
 * on the spacial interval [-5, 10].
 *
 * @param u0 initial data matching std::function<double(double)> and taking
 * in the interval [0,1]
 * @param N number of cells for spacial discretiziation
 * @return vector containing the total masses on an equidistant time grid in
 * [0,3]
 */
/* SAM_LISTING_BEGIN_2 */
template <typename U0Functor>
Eigen::VectorXd traceMass(U0Functor&& u0, unsigned int N) {
  // Spacial boundaries
  double a = -5.0;
  double b = 10.0;

  // Compute \Blue{$\mu(t=0)$}
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, a, b);
  Eigen::VectorXd mu = x.unaryExpr(std::forward<U0Functor>(u0));

#if SOLUTION
  // Get timestep size and number by CLF condition
  double h = (b - a) / N;
  double s_max = std::exp(1.0) - 1.0;  // since \Blue{$0\leq u_0(x)\leq 1$}
  double tau = h / s_max;              // stencil \Blue{$m=1$}
  double T = 3.0;
  // Total number of timesteps
  int M = (int)std::ceil(T / tau);
  tau = 3.0 / M;  // Timestep size

  // Compute solution and total masses at different times
  auto totalMass = [h](const Eigen::VectorXd& mu) { return (mu * h).sum(); };
  Eigen::VectorXd m(M + 1);
  for (int i = 0; i < M; ++i) {
    m(i) = totalMass(mu);
    auto s = [](double u) { return -u; };
    mu = mu + tau * (-1.0 / h) * fluxdiffsource(mu, &godnfn, s, h);
  }
  m(M) = totalMass(mu);

  return m;
#else
  //====================
  // Your code goes here.
  // Replace the dummy return value below by the vector of
  // approximatons of the total mass for the nodes of the
  // temporal mesh.
  return Eigen::VectorXd::Zero(42);
  //====================
#endif
}
/* SAM_LISTING_END_2 */

}  // namespace ConsLawWithSource

#endif  // #define CONSLAWWITHSOURCE
