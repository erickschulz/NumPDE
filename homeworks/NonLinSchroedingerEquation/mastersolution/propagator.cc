/**
 * @file propagator.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 04.05.2020
 * @copyright Developed at ETH Zurich
 */

#include "propagator.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <complex>

namespace NonLinSchroedingerEquation {

// KineticPropagator
/* SAM_LISTING_BEGIN_1 */
KineticPropagator::KineticPropagator(const SparseMatrixXd& A,
                                     const SparseMatrixXcd& M, double tau) {
  // Defeats the rationale of expression templates, but acceptable here, because
  // executed only once in the constructor.
  B_plus_ = M + 0.5 * tau * A.cast<std::complex<double>>();
  SparseMatrixXcd B_minus = M - 0.5 * tau * A.cast<std::complex<double>>();
  // This is the expensive step: LU-factorization of a big sparse matrix.
  // Precomputation is essential
  solver_.compute(B_minus);
}

Eigen::VectorXcd KineticPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  // Cheap elimination steps operating on the LU-factors. Effort is almost O(N)
  // thanks to sophisticated fill-in avoiding techniques employed by the sparse
  // solvers.
  return solver_.solve(B_plus_ * mu);
}
/* SAM_LISTING_END_1 */

// InteractionPropagator
/* SAM_LISTING_BEGIN_2 */
InteractionPropagator::InteractionPropagator(double tau) {
  // We know that the discrete evolution operator for the non-linear part of the
  // method-of-lines ODE boils down to componentwise multiplication of the
  // vector with a single phase shift. Thus, this phase shift can be
  // precomputed, here by means of a lambda function
  phase_multiplier_ = [tau](std::complex<double> z) {
    const std::complex<double> i(0, 1);
    return std::exp(-i * tau * std::norm(z)) * z;
  };
}

Eigen::VectorXcd InteractionPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  // Eigen's way of applying a function to all components of a vector.
  return mu.unaryExpr(phase_multiplier_);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
SplitStepPropagator::SplitStepPropagator(const SparseMatrixXd& A,
                                         const SparseMatrixXcd& M, double tau)
    : kineticPropagator_(A, M, 0.5 * tau), interactionPropagator_(tau) {}

Eigen::VectorXcd SplitStepPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  Eigen::VectorXcd nu(mu.size());
  nu = kineticPropagator_(mu);
  nu = interactionPropagator_(nu);
  nu = kineticPropagator_(nu);
  return nu;
}
/* SAM_LISTING_END_3 */

}  // namespace NonLinSchroedingerEquation
