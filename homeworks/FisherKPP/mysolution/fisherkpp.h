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

    //====================
    // Your code goes here
    //====================
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
    //====================
    // Your code goes here
    //====================
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
