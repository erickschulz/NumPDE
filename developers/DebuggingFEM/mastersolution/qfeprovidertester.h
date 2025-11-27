/**
 * @file qfeprovidertester.h
 * @brief NPDE homework DebuggingFEM code
 * @author Oliver Rietmann
 * @date 03/04/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NPDECODES_DEBUGGINGFEM_QFEPROVIDERTESTER_H_
#define NPDECODES_DEBUGGINGFEM_QFEPROVIDERTESTER_H_

#include <lf/assemble/assemble.h>
#include <lf/assemble/assemble_concepts.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "qfeinterpolator.h"

namespace DebuggingFEM {

/* SAM_LISTING_BEGIN_1 */
template <lf::assemble::EntityMatrixProvider ENTITY_MATRIX_PROVIDER>
class QFEProviderTester {
 public:
  /**
   * @brief Sets up and stores the Galerkin matrix
   * @param dofh dof handler
   * @param element_matrix_provider an instance of LocalLaplaceQFEX, X=1,2,3
   */
  QFEProviderTester(const lf::assemble::DofHandler& dofh,
                    ENTITY_MATRIX_PROVIDER& element_matrix_provider);

  /**
   * @brief Computes the energy (H_1-seminorm) of u
   * @param u function of type double(Eigen::Vector2d)
   */
  template <typename FUNCTOR>
  double energyOfInterpolant(FUNCTOR&& u) const;

 private:
  const lf::assemble::DofHandler& dofh_;
  ENTITY_MATRIX_PROVIDER& element_matrix_provider_;
  Eigen::SparseMatrix<double> A_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <lf::assemble::EntityMatrixProvider ENTITY_MATRIX_PROVIDER>
QFEProviderTester<ENTITY_MATRIX_PROVIDER>::QFEProviderTester(
    const lf::assemble::DofHandler& dofh,
    ENTITY_MATRIX_PROVIDER& element_matrix_provider)
    : dofh_(dofh), element_matrix_provider_(element_matrix_provider) {
#if SOLUTION
  // Assemble the Galerkin matrix and store it into A_
  const lf::base::size_type N_dofs(dofh.NumDofs());
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      mat);
  A_ = Eigen::SparseMatrix<double>(mat.makeSparse());
#else
  //====================
  // Your code goes here
  // Assemble the Galerkin matrix into A_
  //====================
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <lf::assemble::EntityMatrixProvider ENTITY_MATRIX_PROVIDER>
template <typename FUNCTOR>
double QFEProviderTester<ENTITY_MATRIX_PROVIDER>::energyOfInterpolant(
    FUNCTOR&& u) const {
  double energy = 0.0;
#if SOLUTION
  Eigen::VectorXd eta = DebuggingFEM::interpolateOntoQuadFE(dofh_, u);
  // Compute $\cob{\vec{\etabf}^{\top}\VA\vec{\etabf}}$
  energy = eta.dot(A_ * eta);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return energy;
}

/* SAM_LISTING_END_3 */

}  // namespace DebuggingFEM

#endif
