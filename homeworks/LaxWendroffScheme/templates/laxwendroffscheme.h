/**
 * @file laxwendroffscheme.h
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace LaxWendroffScheme {

/**
 * @brief Returns the x-values as needed in the exercise.
 * @param c speed of propagation
 * @param T endtime, T > 0
 * @param M number of time steps, M > 0
 * @return [floor(-3T/h), ..., ceil((3T+1)/h)], where h = tau / e, tau = T / M,
 * e = 2.7182...
 */
Eigen::VectorXd getXValues(double T, unsigned int M);

/**
 * @brief Computes the solution u(x, T) the scalar conservation law with f() =
 * exp() on the nodes x in getXValues(T, M), using the Lax-Wendroff scheme.
 * @param u0 initial data u(x, 0), with x = getXValues(T, M).
 * @param T endtime, T > 0
 * @param M number of time steps, M > 0
 * @return approximation of u(x, T)
 */
Eigen::VectorXd solveLaxWendroff(const Eigen::VectorXd& u0, double T,
                                 unsigned int M);

/**
 * @brief Computes the L1-error of solveLaxWendroff(u0, T, M(i)) against the
 * analytical solution for different values M(i), u0 being the characteristic
 * function of [0, infty) and T = 1.0.
 * @param M vector of numbers of timesteps, M(i) > 0.
 * @return vector of same size as M containing the L1-error for each M(i)
 */
Eigen::VectorXd numexpLaxWendroffRP(const Eigen::VectorXi& M);

/**
 * @brief Provides a reference solution evaluated at the nodes x, obtained by
 * solveLaxWendroff(u0, T, M), where u0 is the smooth inital data from the
 * exercise, T = 1.0 and M = 3200.
 * @param equidistant x-values
 * @return vector of same size as x containing the solution evaluated at x.
 */
Eigen::VectorXd referenceSolution(const Eigen::VectorXd& x);

/**
 * @brief Same as numexpLaxWendroffRP(...), but with the smooth initial data u0
 * from the exercise.
 * @param M vector of numbers of timesteps, M(i) > 0.
 * @return vector of same size as M containing the L1-error for each M(i)
 */
Eigen::VectorXd numexpLaxWendroffSmoothU0(const Eigen::VectorXi& M);

/**
 * @brief Computes the solution u(x, T) the scalar conservation law with f() =
 * exp() on the nodes x in getXValues(T, M), using the Godunov numerical flux
 * (semi-discretization) and explicit euler (full discretization).
 * @param initial data u(x, 0), with x = getXValues(T, M).
 * @param T endtime, T > 0
 * @param M number of time steps, M > 0
 * @return approximation of u(x, T)
 */
Eigen::VectorXd solveGodunov(const Eigen::VectorXd& u0, double T,
                             unsigned int M);

/**
 * @brief Same as numexpLaxWendroffSmoothU0(...), but with solveGodunov(...)
 * instead of solveLaxWendroff(...).
 * @param M vector of numbers of timesteps, M(i) > 0.
 * @return vector of same size as M containing the L1-error for each M(i)
 */
Eigen::VectorXd numexpGodunovSmoothU0(const Eigen::VectorXi& M);

}  // namespace LaxWendroffScheme
