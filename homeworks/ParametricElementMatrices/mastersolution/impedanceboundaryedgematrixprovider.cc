/** @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans, Erick Schulz (refactoring)
 * @date 13/03/2019, 19/11/2019 (refactoring)
 * @copyright Developed at ETH Zurich */

#include "impedanceboundaryedgematrixprovider.h"

namespace ParametricElementMatrices {

/* SAM_LISTING_BEGIN_0 */
ImpedanceBoundaryEdgeMatrixProvider::ImpedanceBoundaryEdgeMatrixProvider(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::VectorXd coeff_expansion)
    : fe_space_(fe_space), coeff_expansion_(coeff_expansion) {
  auto mesh = fe_space->Mesh();
  // Obtain an array of boolean flags for the edges of the mesh: 'true'
  // indicates that the edge lies on the boundary. This predicate will
  // guarantee that the computations are carried only on the boundary edges
  bd_flags_ = std::make_shared<lf::mesh::utils::CodimMeshDataSet<bool>>(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1));
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
bool ImpedanceBoundaryEdgeMatrixProvider::isActive(
    const lf::mesh::Entity& edge) {
  bool is_bd_edge;
  is_bd_edge = (*bd_flags_)(edge);
  return is_bd_edge;
}
/* SAM_LISTING_END_1 */

/**  @brief Compute the local edge element matrix for the Galerkin matrix of
 *
 *           \int_{\boundary \Omega} w(x)^2 u(x) v(x) dx
 *
 * @param edge current edge */
/* SAM_LISTING_BEGIN_2 */
Eigen::MatrixXd ImpedanceBoundaryEdgeMatrixProvider::Eval(
    const lf::mesh::Entity& edge) {
  Eigen::MatrixXd element_matrix(2, 2);

  /* TOOLS AND DATA */
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_->LocGlobMap()};
  // Obtain edge data
  auto edge_global_idx = dofh.GlobalDofIndices(edge);

  // I. COMPUTE LOCAL INTEGRATION DATA
  Eigen::Vector2d w;
  // Evaluate w(x) at the degrees of freedom (endpoints of the edge)
  for (int i = 0; i < 2; i++) {
    w(i) = coeff_expansion_(edge_global_idx[i]);
  }

  // II. COMPUTE LOCAL INTEGRATION DATA
  double edge_length = lf::geometry::Volume(*edge.Geometry());
  // It is assumed here that the edge is straight!
  Eigen::MatrixXd m_1(2, 2), m_2(2, 2), m_3(2, 2);
  m_1 << 24, 6, 6, 4;
  m_2 << 6, 4, 4, 6;
  m_3 << 4, 6, 6, 24;

  // III. PERFORM NUMERICAL QUADRATURE
  element_matrix =
      (w(0) * w(0)) * m_1 + (w(0) * w(1)) * m_2 + (w(1) * w(1)) * m_3;
  element_matrix *= edge_length / 120.;

  return element_matrix;
}
/* SAM_LISTING_END_2 */

}  // namespace ParametricElementMatrices
