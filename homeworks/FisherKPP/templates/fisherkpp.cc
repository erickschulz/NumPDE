/** @file fisherkpp.cc
 *  @brief Homework Problem Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 12.05.20
 *  @copyright Developed at SAM, ETH Zurich
 */

#include "fisherkpp.h"

namespace FisherKPP {

/* Function for the assembly of both Galerkin Matrices,
 * the Mass matrix and the Stiffness matrix.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename DIFF_COEFF>
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler& dofh, DIFF_COEFF&& c) {
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> A_M;
  //====================
  // Your code for matrix assembly goes here
  //====================
  return A_M;
}
/* SAM_LISTING_END_1 */

/* Constructor for StrangSplit */
/* SAM_LISTING_BEGIN_2 */
template <typename DIFF_COEFF>
StrangSplit::StrangSplit(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    double T, unsigned m, double lambda, DIFF_COEFF&& c)
    : fe_space_(fe_space),
      T_(T),
      m_(m),
      lambda_(lambda),
      xi_(1.0 - 0.5 * sqrt(2.0)) {
  const lf::assemble::DofHandler& dofh{fe_space_->LocGlobMap()};
  //====================
  // Your code goes here: initialization of data members
  //====================
}
/* SAM_LISTING_END_2 */

} /* namespace FisherKPP */
