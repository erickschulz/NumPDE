// * Demonstration code for course Numerical Methods for Partial Differential
// * Equations Author: R. Hiptmair, SAM, ETH Zurich
// * Date: May 2020

#ifndef SLFD_H
#define SLFD_H

namespace ConsFV {

// \ref{cpp:fluxdifflc} in lecture document

/* SAM_LISTING_BEGIN_0 */
// arguments:
//   double {h}: meshwidth of equidistant spatial grid
//   Vector {mu}: (finite) vector \Blue{$\vec{\mubf}$} of cell averages
//   Functor {F}: 2-point numerical flux function \Blue{$F=F(v,w)$}
//   Functor {slope}: slope reconstruction function
//   \Blue{$\cor{h\cdot}\sigma_j=\operatorname{slopes}(\mu_{j-1},\mu_j,\mu_{j+1})$}
//
// returns a vector with differences of numerical fluxes,
// which supplies the right-hand side for the FV-MOL ODE
//
// Function that realizes the right hand side operator \Blue{$\Cl_h$} for the
// ODE \eqref{eq:sdev} arising from conservative finite volume
// semidiscretization of the Cauchy problem for a 1D scalar
// conservation law \eqref{eq:clcp}.
template <typename FunctionF, typename FunctionSlopes>
Eigen::VectorXd slopelimfluxdiff(const Eigen::VectorXd& mu, FunctionF&& F,
                                 FunctionSlopes&& slopes) {
  unsigned n = mu.size();  // Number of active dual grid cells
  Eigen::VectorXd sigma = Eigen::VectorXd::Zero(n);  // Vector of slopes
  Eigen::VectorXd fd = Eigen::VectorXd::Zero(n);

  // Computation of slopes \Blue{$\sigma_j$}, uses \Blue{$\mu_0=\mu_1$},
  // \Blue{$m_{N+1}=\mu_N$}, which amounts to constant extension of states
  // beyond domain of influence \Blue{$[a,b]$} of non-constant intial data. Same
  // technique has been applied in \cref{cpp:fluxdiff}
  sigma[0] = slopes(mu[0], mu[0], mu[1]);
  for (unsigned j = 1; j < n - 1; ++j)
    sigma[j] = slopes(mu[j - 1], mu[j], mu[j + 1]);
  sigma[n - 1] = slopes(mu[n - 2], mu[n - 1], mu[n - 1]);

  // Compute linear reconstruction at endpoints of dual cells \eqref{eq:slopval}
  Eigen::VectorXd nup = mu + 0.5 * sigma;
  Eigen::VectorXd num = mu - 0.5 * sigma;

  // As in \cref{cpp:consformevl}: employ constant continuation of data
  // outside \Blue{$[a,b]$}!
  fd[0] = F(nup[0], num[1]) - F(mu[0], num[0]);
  for (unsigned j = 1; j < n - 1; ++j)
    fd[j] =
        F(nup[j], num[j + 1]) - F(nup[j - 1], num[j]);  // see \eqref{eq:2pcf}
  fd[n - 1] = F(nup[n - 1], mu[n - 1]) - F(nup[n - 2], num[n - 1]);
  return fd;
}
/* SAM_LISTING_END_0 */

}  // namespace ConsFV

#endif
