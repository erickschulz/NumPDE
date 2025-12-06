/**
 * @file fvmisentropiceuler.h
 * @brief NPDE homework FVMIsentropicEuler code
 * @author Ralf Hiptmair
 * @date June 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef FVMIsentropicEuler_H_
#define FVMIsentropicEuler_H_

#include <Eigen/Dense>
#include <cstddef>
#include <iostream>

namespace FVMIsentropicEuler {

/** @brief HLLE numerical flux for isentropic Euler equations
 *
 * @tparam functor double -> double
 *
 * @param p pressure function
 * @param pd derivative of pressure function
 * @param v,w the two adjacent cell states
 */
/* SAM_LISTING_BEGIN_1 */
template <typename PFUNCTOR>
Eigen::Vector2d numfluxHLLEIseEul(PFUNCTOR&& p, PFUNCTOR&& pd,
                                  Eigen::Vector2d v, Eigen::Vector2d w) {
  assert((v[0] > 0.0) && "v-density must be positive!");
  assert((w[0] > 0.0) && "w-density must be positive!");
  Eigen::Vector2d nfHLLE{};
  /* **************************************************
     Your code here
     ************************************************** */
  return nfHLLE;
}
/* SAM_LISTING_END_1 */

/** @brief Abstract slope-limited higher order FVM
 *
 * @param mu finite vector of (initial) cell averages
 * @param F functor object for 2-point numerical flux
 * @param slope functor realizing the slope reconstruction function
 * @return right-hand side vector for FVM MOL ODE
 *
 * Function that realizes \$-h\cdot\$ the right hand side operator
 * \$\mathcal{C}_h\$ for the mehtod-of-lines ODE  arising from conservative
 * finite volume semidiscretization of the Cauchy problem for a 1D scalar
 * conservation law.
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FFUNCTOR, typename SLOPEFUNCTOR>
Eigen::VectorXd slopelimfluxdiff_scalar(const Eigen::VectorXd& mu, FFUNCTOR&& F,
                                        SLOPEFUNCTOR&& slopes) {
  size_t n = mu.size();  // Number of active dual grid cells
  Eigen::VectorXd sigma = Eigen::VectorXd::Zero(n);  // Vector of slopes
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);     // Flux differences

  // Computation of slopes \Blue{$\sigma_j$}, uses \Blue{$\mu_0=\mu_1$},
  // \Blue{$\mu_{N+1}=\mu_N$}, which amounts to constant extension of states
  // beyond domain of influence \Blue{$[a,b]$} of non-constant intial data. Same
  // technique has been applied in \lref{cpp:fluxdiff}
  sigma[0] = slopes(mu[0], mu[0], mu[1]);
  for (int j = 1; j < n - 1; ++j)
    sigma[j] = slopes(mu[j - 1], mu[j], mu[j + 1]);
  sigma[n - 1] = slopes(mu[n - 2], mu[n - 1], mu[n - 1]);

  // Compute linear reconstruction at endpoints of dual cells \lref{eq:slopval}
  Eigen::VectorXd nup = mu + 0.5 * sigma;
  Eigen::VectorXd num = mu - 0.5 * sigma;

  // As in \lref{cpp:consformevl}: constant continuation of data outside
  // \Blue{$[a,b]$}!
  fd[0] = F(nup[0], num[1]) - F(mu[0], num[0]);
  for (int j = 1; j < n - 1; ++j)
    // see \lref{eq:2pcf}
    fd[j] = F(nup[j], num[j + 1]) - F(nup[j - 1], num[j]);
  fd[n - 1] = F(nup[n - 1], mu[n - 1]) - F(nup[n - 2], num[n - 1]);
  return fd;
}
/* SAM_LISTING_END_0 */

/** @brief R.h.s. of MOL ODE for abstract high-order slope limited FVM for
 * non-linear systems of conservation laws
 *
 */
/* SAM_LISTING_BEGIN_S */
template <typename FFUNCTOR, typename SLOPEFUNCTOR>
Eigen::MatrixXd slopelimfluxdiff_sys(const Eigen::MatrixXd& mu, FFUNCTOR&& F,
                                     SLOPEFUNCTOR&& slopes) {
  Eigen::Index n = mu.cols();  // Number of active dual grid cells
  // All slopes for all components of cell states
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(mu.rows(), mu.cols());
  Eigen::MatrixXd fd = Eigen::MatrixXd::Zero(mu.rows(), mu.cols());

  /* **************************************************
     Your code here
     ************************************************** */
  return fd;
}
/* SAM_LISTING_END_S */

/** @brief Solve Cauchy problem for 1D isentropic Euler equations by means of
 * minmod-limited MUSCL scheme with HLLE numerical flux
 *
 */
/* SAM_LISTING_BEGIN_F */
template <typename PFUNCTOR, typename RECORDER = std::function<
                                 void(double, const Eigen::MatrixXd&)> >
Eigen::MatrixXd musclIseEul(
    double a, double b, double T, const Eigen::MatrixXd& mu0_mat, PFUNCTOR&& p,
    PFUNCTOR&& dp,
    RECORDER&& rec = [](double /*t*/,
                        const Eigen::MatrixXd& /*mu_mat*/) -> void {}) {
  // Matrix whose columns store the states
  Eigen::MatrixXd mu(mu0_mat.rows(), mu0_mat.cols());
  mu << mu0_mat;
  double h = (b - a) / mu0_mat.cols();
  double t = 0;
  rec(t, mu);
  // Main timestepping loop
  while (t < T) {
    double smax = 0;
    // Adaptive timestep size selection according to (1.2.16)
    for (int j = 0; j < mu.cols(); ++j) {
      double lambda1 = mu(1, j) / mu(0, j) - std::sqrt(dp(mu(0, j)));
      double lambda2 = mu(1, j) / mu(0, j) + std::sqrt(dp(mu(0, j)));
      smax = std::max({smax, std::abs(lambda1), std::abs(lambda2)});
    }
    double tau = std::min(h / (2. * smax), T - t);
    // Lambda function providing HLLE numerical flux
    std::function<Eigen::Vector2d(Eigen::Vector2d, Eigen::Vector2d)> F =
        [p, dp](Eigen::Vector2d u, Eigen::Vector2d v) -> Eigen::Vector2d {
      return numfluxHLLEIseEul(p, dp, u, v);
    };
    // Lambda function for local minmod-limited slope recostruction
    std::function<Eigen::Vector2d(Eigen::Vector2d, Eigen::Vector2d,
                                  Eigen::Vector2d)>
        slopes = [](Eigen::Vector2d u1, Eigen::Vector2d u2,
                    Eigen::Vector2d u3) -> Eigen::Vector2d {
      Eigen::Vector2d out;
      /* **************************************************
         Your code here
         ************************************************** */
      return out;
    };
    // Second-order explicit Heun timestepping based on slopelimfluxdiff\_sys()
    /* **************************************************
       Your code here
       ************************************************** */
    t += tau;
    rec(t, mu);
  }
  return mu;
}
/* SAM_LISTING_END_F */

}  // namespace FVMIsentropicEuler

#endif
