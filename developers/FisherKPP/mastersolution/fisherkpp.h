/** @file fisherkpp.cc
 *  @brief Homework Problem Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 12.05.20
 *  @copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <memory>
#include <utility>

namespace FisherKPP {

/* Function for the assembly of both Galerkin Matrices,
 * the Mass matrix and the Stiffness matrix.
 */
template <typename DIFF_COEFF>
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler& dofh, DIFF_COEFF&& c);

class StrangSplit {
 public:
  // Disabled constructors
  StrangSplit() = delete;
  StrangSplit(const StrangSplit&) = delete;
  StrangSplit(StrangSplit&&) = delete;
  StrangSplit& operator=(const StrangSplit&) = delete;
  StrangSplit& operator=(const StrangSplit&&) = delete;
  // Main constructor
  template <typename DIFF_COEFF>
  explicit StrangSplit(
      const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      double T, unsigned int m, double lambda, DIFF_COEFF&& c);

  // Destructor
  virtual ~StrangSplit() = default;

  /* Member Function StrangSplit
   * Computes the Evolution Operator for the linear parabolic diffusion term
   */
  /* SAM_LISTING_BEGIN_1 */
  Eigen::VectorXd diffusionEvolutionOperator(double tau,
                                             const Eigen::VectorXd& mu) {
    Eigen::VectorXd evol_op;

#if SOLUTION
    // sparse LU decomposition: done in every timestep because timestep size may
    // vary.
    solver.compute(M_ + tau * xi_ * A_);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd rhs = -A_ * mu;

    // First stage SDIRK-2
    Eigen::VectorXd k1 = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    // Second stage SDIRK-2
    Eigen::VectorXd k2 = solver.solve(rhs - tau * (1 - xi_) * A_ * k1);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    // Recover solution
    evol_op = mu + tau * (1 - xi_) * k1 + tau * xi_ * k2;
#else
    //====================
    // Your code goes here
    //====================
#endif
    return evol_op;
  }
  /* SAM_LISTING_END_1 */

  /* Member Function StrangSplit
   * Computes the Evolution for m_ timesteps
   */
  /* SAM_LISTING_BEGIN_2 */
  Eigen::VectorXd Evolution(const Eigen::VectorXd& cap,
                            const Eigen::VectorXd& mu) {
    // Obtain dofhandler
    const lf::assemble::DofHandler& dofh{fe_space_->LocGlobMap()};
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

    Eigen::VectorXd sol(N_dofs);
#if SOLUTION
    Eigen::VectorXd sol_cur(N_dofs);
    Eigen::VectorXd sol_next(N_dofs);

    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N_dofs);

    // step size
    double tau = T_ / m_;

    /* Strang Splitting Method
     * First half time step [0, tau/2]: Diffusion
     */
    sol_cur = diffusionEvolutionOperator(tau / 2., mu);

    // Loop through time steps
    for (int i = 1; i < m_ - 1; i++) {
      // Reaction for next full time step
      sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) *
                                              std::exp(-lambda_ * tau));
      sol_cur = sol_next;
      // Diffusion for next full time step
      sol_next = diffusionEvolutionOperator(tau, sol_cur);
      sol_cur = sol_next;
    }

    // Last half time step: Reaction
    sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) *
                                            std::exp(-lambda_ * tau / 2.));
    sol_cur = sol_next;
    sol = sol_cur;

#else
    //====================
    // Your code goes here
    //====================
#endif
    return sol;
  }
  /* SAM_LISTING_END_2 */

  /* SAM_LISTING_BEGIN_3 */
 private:
  // Finite Element Space
  const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
  // Final Time
  double T_;
  // Number of Timesteps
  unsigned int m_;
  // Growth Factor
  double lambda_;
  // coefficient for SDIRK-2 Butcher Tableau
  double xi_;
  // Galerkin matrix corresponding to the negative Laplacian with Robin Boundary
  // Conditions */
  Eigen::SparseMatrix<double> A_;
  // Galerkin matrix for the Mass Matrix
  Eigen::SparseMatrix<double> M_;
  // Precompute LU decomposition needed for time stepping
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  /* SAM_LISTING_END_3 */
};

} /* namespace FisherKPP. */
