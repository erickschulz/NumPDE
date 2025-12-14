#include <lf/mesh/mesh.h>

#include <Eigen/SparseCore>
#include <memory>

namespace IncidenceMatrices {

/** @brief Create the mesh consisting of a triangle and quadrilateral
 *         from the exercise sheet.
 * @return Shared pointer to the hybrid2d mesh.
 */
std::shared_ptr<lf::mesh::Mesh> createDemoMesh();

/** @brief Compute the edge-vertex incidence matrix G for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
 */
Eigen::SparseMatrix<int> computeEdgeVertexIncidenceMatrix(
    const lf::mesh::Mesh& mesh);

/** @brief Compute the cell-edge incidence matrix D for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
 */
Eigen::SparseMatrix<int> computeCellEdgeIncidenceMatrix(
    const lf::mesh::Mesh& mesh);

/** @brief For a given mesh test if the product of cell-edge and edge-vertex
 *         incidence matrix is zero: D*G == 0?
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return true, if the product is zero and false otherwise
 */
bool testZeroIncidenceMatrixProduct(const lf::mesh::Mesh& mesh);

/** @brief Compute the Hodge-Laplacian matrix L for the given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The Hodge-Laplacian matrix in sparse matrix format
 */
Eigen::SparseMatrix<int> computeHodgeLaplaceMatrix(const lf::mesh::Mesh& mesh);

}  // namespace IncidenceMatrices
