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
  //====================
  // Your code goes here
  //====================
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

  //====================
  // Your code goes here
  //====================

  return element_matrix;
}
/* SAM_LISTING_END_2 */

}  // namespace ParametricElementMatrices
