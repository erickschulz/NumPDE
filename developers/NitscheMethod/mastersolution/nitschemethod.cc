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
#if SOLUTION
  // Fetch geometry object for current cell
  const lf::geometry::Geometry& K_geo{*(cell.Geometry())};
  LF_ASSERT_MSG(K_geo.DimGlobal() == 2, "Mesh must be planar");
  // Obtain physical coordinates of barycenter of triangle
  const Eigen::Vector2d center{K_geo.Global(c_hat_).col(0)};
  // Compute gradients of barycentric coordinate functions
  // Transformation matrix for gradients on reference triangle
  const Eigen::Matrix2d JinvT(K_geo.JacobianInverseGramian(c_hat_));
  // Transform gradients
  const Eigen::Matrix<double, 2, 3> G = JinvT * G_hat_;
  // Element matrix for linear finite elements and the Laplacian
  el_mat = lf::geometry::Volume(K_geo) * G.adjoint() * G;
#else
  //====================
  // Your code goes here
  //====================
#endif
  // ----------------------------------------
  // II: Boundary parts of the bilinear form
  // ----------------------------------------
#if SOLUTION
  // Retrieve pointers to all edges of the triangle
  std::span<const lf::mesh::Entity* const> edges{cell.SubEntities(1)};
  LF_ASSERT_MSG(edges.size() == 3, "Triangle must have three edges!");
  // Loop over edges and check whether they are located on the bondary
  for (int k = 0; k < 3; ++k) {
    if (bd_flags_(*edges[k])) {
      // II(i): Contributions of consistency boundary parts of the bilinear form
      // Edge with local index k is an edge on the boundary
      // Fetch the coordinates of its endpoints
      const lf::geometry::Geometry& ed_geo{*(edges[k]->Geometry())};
      const Eigen::MatrixXd ed_pts{lf::geometry::Corners(ed_geo)};
      // Direction vector of the edge
      const Eigen::Vector2d dir = ed_pts.col(1) - ed_pts.col(0);
      // Rotate counterclockwise by 90 degrees
      const Eigen::Vector2d ed_normal = Eigen::Vector2d(dir(1), -dir(0));
      // For adjusting direction of normal so that it points into the exterior
      // of the domain
      const int ori = (ed_normal.dot(center - ed_pts.col(0)) > 0) ? -1 : 1;
      const Eigen::Matrix3d ed_mat =
          -L_.col(k) * ((ori * ed_normal.transpose()) * G);
      el_mat += (ed_mat + ed_mat.transpose());
      // II(ii): contributions of boundary penalty term
      const double fac = c_ * dir.norm();
      const int l = (k + 1) % 3;
      el_mat(k, k) += fac * 1.0 / 3.0;
      el_mat(k, l) += fac * 1.0 / 6.0;
      el_mat(l, k) += fac * 1.0 / 6.0;
      el_mat(l, l) += fac * 1.0 / 3.0;
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
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
