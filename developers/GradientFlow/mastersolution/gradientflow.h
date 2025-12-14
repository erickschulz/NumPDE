#ifndef GRADIENTFLOW_H_
#define GRADIENTFLOW_H_

/**
 * @file gradientflow.h
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <array>
#include <vector>

namespace GradientFlow {

// Compute the Buther scheme of the SDIRK scheme
Eigen::MatrixXd ButcherMatrix();

/* SAM_LISTING_BEGIN_0 */
// Use Newton method to approximate a stage.
template <typename Functor, typename Jacobian>
Eigen::VectorXd SolveGenStageEquation(Functor&& f, Jacobian&& df,
                                      const Eigen::VectorXd& y,
                                      const Eigen::VectorXd& b, double h,
                                      double rtol = 1E-6, double atol = 1E-8) {
  // Need to solve the equation lhs(g) = g - h*f(y+g)/4 - b = 0.
  // lhs and its Jacobian Jlhs
  auto lhs = [f, y, b, h](const Eigen::VectorXd& g) {
    Eigen::VectorXd val = g - 0.25 * h * f(y + g) - b;
    return val;
  };
  auto Jlhs = [df, y, h](const Eigen::VectorXd& g) {
    // Jlhs(g) = Id - h*df(y+g)/4
    int dim = y.size();
    Eigen::MatrixXd Jval =
        Eigen::MatrixXd::Identity(dim, dim) - 0.25 * h * df(y + g);
    return Jval;
  };

  // Perform Newton iterations:
  Eigen::VectorXd g = Eigen::VectorXd::Zero(y.size());  // initial guess g=0.
  Eigen::VectorXd delta =
      -Jlhs(g).lu().solve(lhs(g));  // Newton correction term.
  int iter = 0;
  int maxiter = 100;  // If correction based termination does not work.
  while (delta.norm() > atol && delta.norm() > rtol * g.norm() &&
         iter < maxiter) {
    g = g + delta;
    delta = -Jlhs(g).lu().solve(lhs(g));
    iter++;
  }
  return g + delta;  // Perform the final step
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
// Compute the stages [g_1, ... , g_5] of the SDIRK method based on Newtons
// method
template <typename Func, typename Jac>
std::array<Eigen::VectorXd, 5> ComputeStages(Func&& f, Jac&& df,
                                             const Eigen::VectorXd& y, double h,
                                             double rtol = 1E-6,
                                             double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> G;  // array of stages
  int dim = y.size();
  Eigen::MatrixXd coeffs = ButcherMatrix();

  for (int i = 0; i < 5; i++) {
    // Calculate coefficients from previous stages:
    Eigen::VectorXd b = Eigen::VectorXd::Zero(dim);
    for (int j = 0; j < i; j++) {
      b += coeffs(i, j) * f(y + G[j]);
    }
    b *= h;
    // Compute stage i and store in G.at(i):
    G[i] = SolveGenStageEquation(f, df, y, b, h, rtol, atol);
  }

  return G;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
// Compute one step of the SDIRK scheme
template <typename Func, typename Jac>
Eigen::VectorXd DiscEvolSDIRK(Func&& f, Jac&& df, const Eigen::VectorXd& y,
                              double h, double rtol = 1E-6,
                              double atol = 1E-8) {
  // The b weights are in the last row of coeffs.
  Eigen::MatrixXd coeffs = ButcherMatrix();
  int n_stages = coeffs.cols();
  Eigen::VectorXd b = coeffs.row(n_stages);

  // Vector after one step of the SDIRK scheme
  Eigen::VectorXd Psi;

#if SOLUTION
  Psi = y;
  // Compute array of stages
  std::array<Eigen::VectorXd, 5> G = ComputeStages(f, df, y, h, rtol, atol);
  for (int i = 0; i < n_stages; i++) {
    Psi += h * b(i) * f(y + G[i]);
  }

#else
  //====================
  // Your code goes here
  //====================
#endif
  return Psi;
}
/* SAM_LISTING_END_4 */
// Solve the gradient flow problem based on the SDIRK scheme using M uniform
// timesteps. Return the full approximated solution trajectory
std::vector<Eigen::VectorXd> SolveGradientFlow(const Eigen::VectorXd& d,
                                               double lambda,
                                               const Eigen::VectorXd& y0,
                                               double T, unsigned int M);

}  // namespace GradientFlow

#endif  // #ifndef GRADIENTFLOW_H_
