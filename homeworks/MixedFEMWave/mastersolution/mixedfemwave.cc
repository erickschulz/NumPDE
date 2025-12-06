/**
 * @file mixedfemwave.cc
 * @brief NPDE homework MixedFEMWave
 * @author Erick Schulz
 * @date 14.06.2020
 * @copyright Developed at ETH Zurich
 */

#include "mixedfemwave.h"

namespace MixedFEMWave {

/* SAM_LISTING_BEGIN_A */
// An auxiliary function providing a QuadRule object for the local trapezoidal
// rule on a triangle
lf::quad::QuadRule make_TriaQR_TrapezoidalRule() {
  // A 3-point quadrature rule
  Eigen::MatrixXd points(2, 3);
  Eigen::VectorXd weights(3);
  // Quadrature nodes, reference coordinates of vertices of reference triangle
  points(0, 0) = 0.0;
  points(1, 0) = 0.0;
  points(0, 1) = 1.0;
  points(1, 1) = 0.0;
  points(0, 2) = 0.0;
  points(1, 2) = 1.0;
  // Quadrature weights on reference triangle
  weights(0) = 0.16666666666666665741;
  weights(1) = 0.16666666666666665741;
  weights(2) = 0.16666666666666665741;
  // Quadrature rule is exact for polynomials of degree 1
  return lf::quad::QuadRule(lf::base::RefEl::kTria(), std::move(points),
                            std::move(weights), 1);
}
/* SAM_LISTING_END_A */

/* SAM_LISTING_BEGIN_5 */
// Auxiliary function: Determine combined areas of cells adjacent to the nodes
// of a mesh
lf::mesh::utils::CodimMeshDataSet<double> areasOfAdjacentCells(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  const lf::mesh::Mesh& mesh{*mesh_p};
  // The areas are stored in a node-indexed (codim == 2) array
  lf::mesh::utils::CodimMeshDataSet<double> areas(mesh_p, 2, 0.0);
  // Loop over all cells (codim == 0 entities)
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    const double area = lf::geometry::Volume(*(cell->Geometry()));
    // Loop over nodes of  a cell (relative codim == 2 sub-entities)
    const std::span<const lf::mesh::Entity* const> cell_nodes{
        cell->SubEntities(2)};
    for (const lf::mesh::Entity* node : cell_nodes) {
      LF_ASSERT_MSG(node->RefEl() == lf::base::RefEl::kPoint(),
                    "Illegal topological type for node");
      areas(*node) += area;
    }
  }
  return areas;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double> computeMQ(const lf::assemble::DofHandler& dofh_Q) {
  // TOOLS AND DATA
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs_Q(dofh_Q.NumDofs());
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh_Q.Mesh();
  // Sparse matrix to be returned
  Eigen::SparseMatrix<double> M_Q(N_dofs_Q, N_dofs_Q);
  // ASSEMBLY
  // First option: assembly via triplet format
  lf::assemble::COOMatrix<double> M_Q_COO(N_dofs_Q, N_dofs_Q);
  // Loop over all cells
  for (const lf::mesh::Entity* entity : mesh_p->Entities(0)) {
    const double area = lf::geometry::Volume(*(entity->Geometry()));
    auto global_idx = dofh_Q.GlobalDofIndices(*entity);
    M_Q_COO.AddToEntry(global_idx[0], global_idx[0], area);
    M_Q_COO.AddToEntry(global_idx[1], global_idx[1], area);
  }
  M_Q = M_Q_COO.makeSparse();
  std::cout << "Assembly: M_Q finished" << std::endl;
  // Allternative option: reserve() method plus direct initialization of
  // diagonal entries.

  return M_Q;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_9 */
Eigen::SparseMatrix<double> computeB(const lf::assemble::DofHandler& dofh_V,
                                     const lf::assemble::DofHandler& dofh_Q) {
  // TOOLS AND DATA
  auto mesh_p = dofh_V.Mesh();
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs_Q(dofh_Q.NumDofs());
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());
  // Sparse matrix to be returned
  Eigen::SparseMatrix<double> B(N_dofs_Q, N_dofs_V);

  // ASSEMBLY
  lf::assemble::COOMatrix<double> B_COO(N_dofs_Q, N_dofs_V);
  BElemMatProvider BElemMat_builder;
  lf::assemble::AssembleMatrixLocally(0, dofh_V, dofh_Q, BElemMat_builder,
                                      B_COO);

  // Node-index array of flags marking boundary nodes
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Set columns of the matrix to zero using a predicate that evaluates to true,
  // if the column index is that of a d.o.f. associated with a node on the
  // boundary
  B_COO.setZero([&bd_flags, &dofh_V](lf::assemble::gdof_idx_t /*i*/,
                                     lf::assemble::gdof_idx_t j) -> bool {
    return bd_flags(dofh_V.Entity(j));
  });

  // Convert from triplet to CRS format
  B = B_COO.makeSparse();

  return B;
}
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_B */
Eigen::Matrix<double, 2, 3> BElemMatProvider::Eval(
    const lf::mesh::Entity& tria) {
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const auto vertices = lf::geometry::Corners(*(tria.Geometry()));
  LF_ASSERT_MSG((vertices.cols() == 3) && (vertices.rows() == 2),
                "Invalid vertex coordinate " << vertices.rows() << "x"
                                             << vertices.cols() << " matrix");
  // Matrix for returning result
  Eigen::Matrix<double, 2, 3> elmat;
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix grad\_bary\_coords
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();  // area of triangular cell
  const double area = 0.5 * std::abs(X.determinant());
  // This matrix contains the barycentric coordinate functions in its columns
  auto inv = X.inverse().eval();
  const Eigen::Matrix<double, 2, 3> grad_bary_coords{inv.block<2, 3>(1, 0)};
  // Since the local shape function for the finite element space $Q_h$ are
  // Cartesian coordinate vectors, we just need to scale the components of the
  // gradients of the barycentric coordinate functions with the area of the
  // triangle.
  elmat = area * grad_bary_coords;
  return elmat;
}
/* SAM_LISTING_END_B */

}  // namespace MixedFEMWave
