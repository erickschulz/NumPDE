#ifndef SLOPELIMFLUXDIFF_H_
#define SLOPELIMFLUXDIFF_H_

/**
 * @file slopelimfluxdiff.h
 * @brief NPDE homework ExtendedMUSCL code
 * @author Ralf Hiptmair
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace ExtendedMUSCL {

/* SAM_LISTING_BEGIN_0 */
// Function parameters:
//   Vector texttt{mu}: (finite) vector \Blue{$\vec{\mubf}$} of cell averages
//   Functor texttt{F}: 2-point numerical flux function \Blue{$F=F(v,w)$}
//   Functor texttt{slope}: slope reconstruction function
//   \Blue{$h\sigma_j=\operatorname{slopes}(\mu_{j-1},\mu_j,\mu_{j+1})$}
//
// returns a vector with differences of numerical fluxes
//
// Function that realizes \cob{$-h\cdot$} the right hand side operator
// \Blue{$\Cl_h$} for the mehtod-of-lines ODE  arising from conservative finite
// volume semidiscretization of the Cauchy problem for a 1D scalar conservation
// law \lref{eq:clcp}.
template <typename FunctionF, typename FunctionSlopes>
Eigen::VectorXd slopelimfluxdiff(const Eigen::VectorXd& mu, FunctionF&& F,
                                 FunctionSlopes&& slopes) {
  int n = mu.size();  // Number of active dual grid cells
  Eigen::VectorXd sigma = Eigen::VectorXd::Zero(n);  // Vector of slopes
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);     // Flux differences

  // Computation of slopes \Blue{$\sigma_j$}, uses \Blue{$\mu_0=\mu_1$},
  // \Blue{$m_{N+1}=\mu_N$}, which amounts to constant extension of states
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

/* SAM_LISTING_BEGIN_1 */
// Function parameters:
//   double h: meshwidth of equidistant spatial grid
//   Vector mu: (finite) vector \Blue{$\vec{\mubf}$} of cell averages
//   Functor F: 2-point numerical flux function \Blue{$F=F(v,w)$}
//   Functor slope: slope reconstruction function
//   \Blue{$h\sigma_j=\operatorname{slopes}(\mu_{j-1},\mu_j,\mu_{j+1})$}
//
// returns a vector with differences of numerical fluxes
//
// Function that realizes \cob{$-h\cdot$} the right hand side operator
// \Blue{$\Cl_h$} for the mehtod-of-lines ODE  arising from conservative finite
// volume semidiscretization of the Cauchy problem for a 1D scalar conservation
// law \lref{eq:clcp} in a \samemp*{1-periodic setting}.
template <typename FunctionF, typename FunctionSlopes>
Eigen::VectorXd slopelimfluxdiffper(const Eigen::VectorXd& mu, FunctionF&& F,
                                    FunctionSlopes&& slopes) {
  int n = mu.size();  // Number of active dual grid cells
  Eigen::VectorXd sigma = Eigen::VectorXd::Zero(n);  // Vector of slopes
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);     // Flux differences

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
/* SAM_LISTING_END_1 */

}  // namespace ExtendedMUSCL

#endif  // SLOPELIMFLUXDIFF_H_
