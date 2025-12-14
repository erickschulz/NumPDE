/**
 * @ file
 * @ brief NPDE homework on Nitsche's method
 * @ author R. Hiptmair
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace NitscheMethod {

/**
 * @brief Flawed edge based matrix provider for boundary contributions
 *
 *  This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 */
/* SAM_LISTING_BEGIN_2 */
class NitscheBoundaryMatProvider {
 public:
  NitscheBoundaryMatProvider(lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags,
                             double c)
      : bd_flags_(bd_flags), c_(c) {}
  virtual ~NitscheBoundaryMatProvider() = default;
  [[nodiscard]] virtual bool isActive(const lf::mesh::Entity& edge) const {
    return bd_flags_(edge);
  }
  [[nodiscard]] Eigen::Matrix2d Eval(const lf::mesh::Entity& edge) const;

 private:
  lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags_;
  double c_;
};
/* SAM_LISTING_END_2 */

/**
 * @brief Element matrix builder class for Nitsche's method
 *
 *  This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 */
/* SAM_LISTING_BEGIN_3 */
class LinearFENitscheElementMatrix {
 public:
  LinearFENitscheElementMatrix(
      lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags, double c)
      : bd_flags_(bd_flags), c_(c) {}
  virtual ~LinearFENitscheElementMatrix() = default;
  [[nodiscard]] virtual bool isActive(const lf::mesh::Entity& /*cell*/) const {
    return true;
  }
  [[nodiscard]] Eigen::Matrix3d Eval(const lf::mesh::Entity& cell) const;

 private:
  lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags_;
  double c_;
  // Constant matrices for use in Eval()
  const Eigen::MatrixXd c_hat_ = Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0);
  const Eigen::MatrixXd G_hat_ =
      (Eigen::Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0)
          .finished();
  const Eigen::Matrix3d L_ =
      (Eigen::Matrix3d() << 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5)
          .finished();
};
/* SAM_LISTING_END_3 */

/**
 * @brief ENTITY_VECTOR_PROVIDER for Nitsche's method
 *
 * @tparam FUNCTOR type equipped with double operator(Vector2d) const
 *
 * Class for computation of element vectors for linear finite elements
 */
/* SAM_LISTING_BEGIN_6 */
template <typename FUNCTOR>
class NitscheElemVecProvider {
 public:
  NitscheElemVecProvider(FUNCTOR g,
                         lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags,
                         double c)
      : g_(g), bd_flags_(bd_flags), c_(c) {}
  virtual ~NitscheElemVecProvider() = default;
  [[nodiscard]] virtual bool isActive(const lf::mesh::Entity& /*cell*/) {
    return true;
  }
  [[nodiscard]] Eigen::Vector3d Eval(const lf::mesh::Entity& cell) const;

 private:
  FUNCTOR g_;  // $\cob{\Bx\mapsto g(\Bx)}$
  // Flags marking edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool>& bd_flags_;
  double c_;  // Penalty parameter $\cob{c}$
  // Reference coordinate of barycentre of a triangle
  const Eigen::MatrixXd c_hat_ = Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0);
  // Gradients of reference shape functions
  const Eigen::MatrixXd G_hat_ =
      (Eigen::Matrix<double, 2, 3>() << -1.0, 1.0, 0.0, -1.0, 0.0, 1.0)
          .finished();
};
/* SAM_LISTING_END_6 */

// Implementation of local computations yielding element vector
/* SAM_LISTING_BEGIN_7 */
template <typename FUNCTOR>
Eigen::Vector3d NitscheElemVecProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& cell) const {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  // Fetch geometry object for current cell
  const lf::geometry::Geometry& K_geo{*(cell.Geometry())};
  const Eigen::MatrixXd cell_pts{lf::geometry::Corners(K_geo)};
  LF_ASSERT_MSG(K_geo.DimGlobal() == 2, "Mesh must be planar");
  // Obtain physical coordinates of barycenter of triangle
  const Eigen::Vector2d center{K_geo.Global(c_hat_).col(0)};
  // Compute gradients of barycentric coordinate functions
  // Transformation matrix for gradients on reference triangle
  const Eigen::Matrix2d JinvT(K_geo.JacobianInverseGramian(c_hat_));
  // Transform gradients
  const Eigen::Matrix<double, 2, 3> G = JinvT * G_hat_;
  // Element vector
  Eigen::Vector3d el_vec = Eigen::Vector3d::Zero();
  // Loop over edges and check whether they are
  // located on the bondary
  std::span<const lf::mesh::Entity* const> edges{cell.SubEntities(1)};
  for (int k = 0; k < 3; ++k) {
    if (bd_flags_(*edges[k])) {
      // Edge with local index k is an edge on the boundary
      // Its one endpoints has local vertex index k
      // Fetch the local index of the other end node
      const int l = (k + 1) % 3;
      // I: Contribution from consistency correction
      // Direction vector of the edge
      const Eigen::Vector2d dir = cell_pts.col(k) - cell_pts.col(l);
      // Rotate counterclockwise by 90 degrees
      const Eigen::Vector2d ed_normal = Eigen::Vector2d(dir(1), -dir(0));
      // For adjusting direction of normal
      const int ori = (ed_normal.dot(center - cell_pts.col(k)) > 0) ? -1 : 1;
      // Evaluate Dirichlet data function in the endpoints
      const double g0 = g_(cell_pts.col(k));
      const double g1 = g_(cell_pts.col(l));
      // Dirichlet data approximated as piecewise linear
      el_vec -= 0.5 * (g0 + g1) * (G.transpose() * (ori * ed_normal));
      // II: Contribution from penalty correction
      const double fac = c_ * dir.norm();
      el_vec[k] += fac * (g0 / 3.0 + g1 / 6.0);
      el_vec[l] += fac * (g1 / 3.0 + g0 / 6.0);
    }
  }
  return el_vec;
}
/* SAM_LISTING_END_7 */

/**
 * @brief Spurious assembly of Galerkin matrix for Nitsches method
 */
Eigen::SparseMatrix<double> assembleNitscheGalerkinMatrix(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> lin_fes_p,
    double c);

/**
 * @brief Correct assembly of Galerkin matrix for Nitsches method
 */
Eigen::SparseMatrix<double> computeNitscheGalerkinMatrix(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> lin_fes_p,
    double c);

/**
 * @brief Assembly of right-hand-side vector for Nitsche's method
 */
template <typename FUNCTOR>
Eigen::VectorXd computeNitscheLoadVector(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> lin_fes_p,
    FUNCTOR g, double c) {
  // Pointer to underlying mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{lin_fes_p->Mesh()};
  // Obtain local-to-global index mapper "D.o.f. handler"
  const lf::assemble::DofHandler& dofh{lin_fes_p->LocGlobMap()};
  // Flag all edge (co-dimension-1 entities) on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Provider object for element vectors
  NitscheElemVecProvider<FUNCTOR> nitsche_vec_builder(g, bd_flags, c);
  // Object for sparse matrix to be filled by cell-oriented assembly
  const int N_dofs = dofh.NumDofs();
  Eigen::VectorXd phi{Eigen::VectorXd::Zero(N_dofs)};
  // Cell-oriented assembly
  lf::assemble::AssembleVectorLocally(0, dofh, nitsche_vec_builder, phi);
  return phi;
}

}  // namespace NitscheMethod
