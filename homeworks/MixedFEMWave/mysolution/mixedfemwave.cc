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
  // =============================================
  // Your code here
  // =============================================

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

  // ========================================
  // Your code here
  // ========================================

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
  // ========================================
  // Your code here
  // ========================================
  return elmat;
}
/* SAM_LISTING_END_B */

}  // namespace MixedFEMWave
