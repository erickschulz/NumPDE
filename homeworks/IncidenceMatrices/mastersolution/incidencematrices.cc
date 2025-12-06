#include "incidencematrices.h"

#include <Eigen/src/SparseCore/SparseUtil.h>
#include <lf/base/base.h>
#include <lf/base/lf_assert.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace IncidenceMatrices {

/** @brief Create the mesh consisting of a triangle and quadrilateral
 *         from the exercise sheet.
 * @return Shared pointer to the hybrid2d mesh.
 */
std::shared_ptr<lf::mesh::Mesh> createDemoMesh() {
  // builder for a hybrid mesh in a world of dimension 2
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Add points
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 0});    // (0)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 0});    // (1)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 1});    // (2)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 1});    // (3)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  Eigen::MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      std::array<lf::mesh::Mesh::size_type, 3>{
          {1, 2, 4}},  // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  Eigen::MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      std::array<lf::mesh::Mesh::size_type, 4>{{0, 1, 4, 3}},
      std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<lf::mesh::Mesh> demoMesh_p = mesh_factory_ptr->Build();

  return demoMesh_p;
}

/** @brief Compute the edge-vertex incidence matrix G for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<int> computeEdgeVertexIncidenceMatrix(
    const lf::mesh::Mesh& mesh) {
  // Store edge-vertex incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> G;

  // Mesh::NumEntities(unsigned codim) returns the number of elements
  // with given codimension. Codim(Edge) = 1, Codim(Node) = 2.
  const lf::mesh::Mesh::size_type numEdges = mesh.NumEntities(1),
                                  numNodes = mesh.NumEntities(2);
  // Following the demo for the reserve()-initialising the sparse matrix given
  // in the exercise sheet. From (2.1a) we know that G has exactly 2 entries
  // per row.
  G = Eigen::SparseMatrix<int, Eigen::RowMajor>(numEdges, numNodes);
  G.reserve(Eigen::VectorXi::Constant(numEdges, 2));

  // To compute G efficiently we iterate over all edges and check the index
  // of the nodes at its end. This is the efficient way to do the assembly,
  // introduced as "distribute scheme" in class. We cannot iterative over
  // vertices, because LehrFEM++ does not allow to visit the edges
  // adjacent to a vertex
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    // Get index of this edge
    lf::mesh::Mesh::size_type edgeIdx = mesh.Index(*edge);
    // Get the nodes and their indices.
    // Note, that seen from the edges the nodes have codim 1, not 2,
    // hence we call SubEntities(1). This is a relative codimension!
    auto nodes = edge->SubEntities(1);
    lf::mesh::Mesh::size_type firstNodeIdx = mesh.Index(*nodes[0]);
    lf::mesh::Mesh::size_type lastNodeIdx = mesh.Index(*nodes[1]);
    // Add the matrix entries according to the definition
    G.coeffRef(edgeIdx, firstNodeIdx) += 1;
    G.coeffRef(edgeIdx, lastNodeIdx) -= 1;
  }

  return G;
}
/* SAM_LISTING_END_1 */

/** @brief Compute the cell-edge incidence matrix D for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<int> computeCellEdgeIncidenceMatrix(
    const lf::mesh::Mesh& mesh) {
  // Store cell-edge incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> D;

  // Mesh::NumEntities(unsigned codim) returns the number of elements
  // with given codimension. Codim(Edge) = 0, Codim(Node) = 1.
  const lf::mesh::Mesh::size_type numCells = mesh.NumEntities(0),
                                  numEdges = mesh.NumEntities(1);
  // Following the demo for the reserve()-initialising the sparse matrix given
  // in the exercise sheet. From (2.1a) we know that D has at most 4 entries
  // per row.
  D = Eigen::SparseMatrix<int, Eigen::RowMajor>(numCells, numEdges);
  D.reserve(Eigen::VectorXi::Constant(numCells, 4));

  // To compute D efficiently we iterate over all cells and check the
  // orientations (+1 or -1, same as in the definition of the matrix D)
  // of its edges. For this we may use RelativeOrientations().
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Get cell index
    lf::mesh::Mesh::size_type cellIdx = mesh.Index(*cell);
    // Get edges and their orientations (these already the entries for D!)
    auto edges = cell->SubEntities(1);
    auto edgeOrientations = cell->RelativeOrientations();

    // Iterate over both and add to D
    auto edgeIt = edges.begin();
    auto orntIt = edgeOrientations.begin();
    for (; edgeIt != edges.end() && orntIt != edgeOrientations.end();
         ++edgeIt, ++orntIt) {
      // Get the edge index and add its orientation to D
      lf::mesh::Mesh::size_type edgeIdx = mesh.Index(**edgeIt);
      D.coeffRef(cellIdx, edgeIdx) += lf::mesh::to_sign(*orntIt);
    }
  }

  return D;
}
/* SAM_LISTING_END_2 */

/** @brief For a given mesh test if the product of cell-edge and edge-vertex
 *        incidence matrix is zero: D*G == 0?
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *             such as lf::mesh::hybrid2d::Mesh)
 * @return true, if the product is zero and false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testZeroIncidenceMatrixProduct(const lf::mesh::Mesh& mesh) {
  bool isZero = false;

  Eigen::SparseMatrix<int> G = computeEdgeVertexIncidenceMatrix(mesh),
                           D = computeCellEdgeIncidenceMatrix(mesh);

  Eigen::SparseMatrix<int> O = D * G;
  // Possibility 1:
  // Not prone to roundoff errors, since O is an integer matrix!
  isZero = O.norm() == 0;
  // Possibility 2: But this doesn't use the fact that O is sparse.
  // isZero = Eigen::MatrixXi(O).isZero(0);
  return isZero;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<int> computeHodgeLaplaceMatrix(const lf::mesh::Mesh& mesh) {
  // Size of Hodge Laplacian matrix for discrete 1-forms is equal to the number
  // of edges of the mesh
  const size_t N = mesh.NumEntities(1);
  // Store cell-edge incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> L(N, N);
  // Triplet vector to be used for the initialization of the sparse matrix
  std::vector<Eigen::Triplet<int>> triplets;
  // Pass I: Cell-oriented assembly of $\VD_1^{\top}\VD_1$.
  // Visit all cells and compute contributions of pairs of faces
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Array of pointers to the edges of the cell, see \lref{ex:subent}
    std::span<const lf::mesh::Entity* const> edges{cell->SubEntities(1)};
    const unsigned int n_ed = edges.size();
    LF_ASSERT_MSG(n_ed <= 4, "Illegal number of edges");
    // Inquire about relative orientations, see \lref{rem:ori}
    std::span<const lf::mesh::Orientation> ror{cell->RelativeOrientations()};
    LF_ASSERT_MSG(ror.size() == n_ed, "Size mismatch for orientation array");
    // Obtain global indices of edges
    std::array<double, 4> edge_idx;
    for (unsigned int j = 0; j < n_ed; ++j) {
      edge_idx[j] = mesh.Index(*edges[j]);
    }
    // Initialize triplets corresponding to pairs of edges
    for (unsigned int j = 0; j < n_ed; ++j) {
      for (unsigned int i = 0; i < n_ed; ++i) {
        triplets.emplace_back(edge_idx[j], edge_idx[i],
                              lf::mesh::to_sign(ror[j]) * to_sign(ror[i]));
      }
    }
  }
  // Pass II: Node-oriented assembly of $\VD_0 \VD_0^{\top}$
  // Preparatory step: build lists of edges adjacent to nodes
  std::vector<std::vector<std::pair<size_t, int>>> eds_node(
      mesh.NumEntities(2));
  // Run through the edges and collect their endpoints
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    // Get index of this edge
    lf::mesh::Mesh::size_type edgeIdx = mesh.Index(*edge);
    // Get the nodes and their indices.
    // Note, that seen from the edges the nodes have codim 1, not 2,
    // hence we call SubEntities(1). This is a relative codimension!
    auto nodes = edge->SubEntities(1);
    lf::mesh::Mesh::size_type firstNodeIdx = mesh.Index(*nodes[0]);
    lf::mesh::Mesh::size_type lastNodeIdx = mesh.Index(*nodes[1]);
    LF_ASSERT_MSG(firstNodeIdx < mesh.NumEntities(2),
                  "First node index out of range");
    LF_ASSERT_MSG(lastNodeIdx < mesh.NumEntities(2),
                  "Last node index out of range");
    eds_node[firstNodeIdx].emplace_back(edgeIdx, -1);
    eds_node[lastNodeIdx].emplace_back(edgeIdx, 1);
  }
  // Run through nodes and pairs ofv adjacent edges
  for (const std::vector<std::pair<size_t, int>>& edgesinfo : eds_node) {
    for (const std::pair<size_t, int>& idxori1 : edgesinfo) {
      for (const std::pair<size_t, int>& idxori2 : edgesinfo) {
        triplets.emplace_back(idxori1.first, idxori2.first,
                              idxori1.second * idxori2.second);
      }
    }
  }
  L.setFromTriplets(triplets.begin(), triplets.end());
  return L;
}
/* SAM_LISTING_END_4 */

}  // namespace IncidenceMatrices
