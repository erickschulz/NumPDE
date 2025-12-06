/**
 * @ file
 * @ brief NPDE homework on Nitsche's method (designed an exam problem)
 * @ author R. Hiptmair
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "nitschemethod.h"

namespace NitscheMethod {

// Do-nothing implementation
Eigen::Matrix2d NitscheBoundaryMatProvider::Eval(
    const lf::mesh::Entity& edge) const {
  LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                "Entity must be an edge");
  // No implementation given here
  return Eigen::Matrix2d::Zero();
}

// Implementation of local computation of element matrix for Nitsche's method
/* SAM_LISTING_BEGIN_5 */
Eigen::Matrix3d LinearFENitscheElementMatrix::Eval(
    const lf::mesh::Entity& cell) const {
  LF_ASSERT_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Only implemented for triangles");
  // The element matrix returned from this function
  Eigen::Matrix3d el_mat = Eigen::Matrix3d::Zero();
  // ------------------------------------------------------------------
  // I: Compute element matrix induced by volume part of bilinear form
  // ------------------------------------------------------------------
  //====================
  // Your code goes here
  //====================
  // ----------------------------------------
  // II: Boundary parts of the bilinear form
  // ----------------------------------------
  //====================
  // Your code goes here
  //====================
  return el_mat;
}  // end Eval()
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_7 */
Eigen::SparseMatrix<double> assembleNitscheGalerkinMatrix(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> lin_fes_p,
    double c) {
  // (Pointer to) underlying mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{lin_fes_p->Mesh()};
  // Obtain local-to-global index mapper
  const lf::assemble::DofHandler& dofh{lin_fes_p->LocGlobMap()};
  // Provider object for element matrices for negative Laplacian and linear
  // Lagrangian finite elements
  lf::uscalfe::LinearFELaplaceElementMatrix lapl_elmat_builder{};
  // Object for sparse matrix to be filled by cell-oriented assembly
  const int N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Invoke cell-oriented assembly for the volume part of the finite-element
  // Galerkin matrix
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, lapl_elmat_builder, A);
  // Flag all edge (co-dimension-1 entities) on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  NitscheBoundaryMatProvider nitsche_bd_builder(bd_flags, c);
  // Add boundary contributions to Galerkin matrix stored in A
  // Assembly covers edges (co-dimensions-1 entities)!
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, nitsche_bd_builder, A);
  // Convert the matrix from triplet format to CRS format
  return A.makeSparse();
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
Eigen::SparseMatrix<double> computeNitscheGalerkinMatrix(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> lin_fes_p,
    double c) {
  // Pointer to underlying mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{lin_fes_p->Mesh()};
  // Obtain local-to-global index mapper "D.o.f. handler"
  const lf::assemble::DofHandler& dofh{lin_fes_p->LocGlobMap()};
  // Flag all edge (co-dimension-1 entities) on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Provider object for element matrices
  LinearFENitscheElementMatrix nitsche_mat_builder(bd_flags, c);
  // Object for sparse matrix to be filled by cell-oriented assembly
  const int N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Cell-oriented assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, nitsche_mat_builder, A);
  return A.makeSparse();
}
/* SAM_LISTING_END_8 */

}  // namespace NitscheMethod
