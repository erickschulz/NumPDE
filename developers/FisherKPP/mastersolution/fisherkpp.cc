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
#if SOLUTION
  // Obtain mesh and finite element space
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofh.Mesh();
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // Matrix format for assembly
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);

  // Function handles for Stiffness and Mass matrix
  auto one = [](Eigen::Vector2d x) -> double { return 1.0; };
  auto zero = [](Eigen::Vector2d x) -> double { return 0.0; };

  // Mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_c{c};
  lf::mesh::utils::MeshFunctionGlobal mf_one{one};
  lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};

  // Element matrix provider
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_c),
                                                      decltype(mf_zero)>
      elMat_Stiff(fe_space, mf_c, mf_zero);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_zero),
                                                      decltype(mf_one)>
      elMat_Mass(fe_space, mf_zero, mf_one);

  // Assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Stiff, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Mass, M_COO);

  // Sparse matrix format
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();

  A_M = std::make_pair(A, M);

#else
  //====================
  // Your code for matrix assembly goes here
  //====================
#endif
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
#if SOLUTION
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      galerkinpair = assembleGalerkinMatrices(dofh, c);
  A_ = galerkinpair.first;
  M_ = galerkinpair.second;
  // Butcher Tableau Coefficient for SDIRK-2

#else
  //====================
  // Your code goes here: initialization of data members
  //====================
#endif
}
/* SAM_LISTING_END_2 */

} /* namespace FisherKPP */
