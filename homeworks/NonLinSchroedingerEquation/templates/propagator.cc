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
  //====================
  // Your code goes here
  //====================
}

Eigen::VectorXcd KineticPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  //====================
  // Your code goes here
  // Replace mu by its value after a timestep tau
  return mu;
  //====================
}
/* SAM_LISTING_END_1 */

// InteractionPropagator
/* SAM_LISTING_BEGIN_2 */
InteractionPropagator::InteractionPropagator(double tau) {
  //====================
  // Your code goes here
  //====================
}

Eigen::VectorXcd InteractionPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  //====================
  // Your code goes here
  // Replace mu by its value after a timestep tau
  return mu;
  //====================
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
//====================
// Your code goes here
// Change this dummy implementation of the constructor:
SplitStepPropagator::SplitStepPropagator(const SparseMatrixXd& A,
                                         const SparseMatrixXcd& M, double tau) {
}
//====================

Eigen::VectorXcd SplitStepPropagator::operator()(
    const Eigen::VectorXcd& mu) const {
  Eigen::VectorXcd nu(mu.size());
  //====================
  // Your code goes here
  // Implement the Strang splitting
  //====================
  return nu;
}
/* SAM_LISTING_END_3 */

}  // namespace NonLinSchroedingerEquation
