/** @file
 * @brief NPDE BoundaryWave
 * @author Erick Schulz
 * @date 24/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "boundarywave.h"

namespace BoundaryWave {

/* SAM_LISTING_BEGIN_1 */
lf::assemble::COOMatrix<double> buildM(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p) {
  // I. TOOLS AND DATA
  // Pointer to current fe_space and mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p(fe_space_p->Mesh());
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // II : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
#if SOLUTION
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Creating a predicate that will guarantee that the computations are carried
  // only on the edges of the mesh using the boundary flags
  // Actually a redundant step, because 'bdflags' is a predicate already.
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge);
  };
  // Coefficient function used in the class template MassEdgeMatrixProvider
  auto eta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(eta),
                                      decltype(edges_predicate)>
      edgemat_builder(fe_space_p, eta, edges_predicate);
  // Invoke assembly on edges by specifying co-dimension = 1
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, M);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return M;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
lf::assemble::COOMatrix<double> buildA(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p) {
  // I. TOOLS AND DATA
  // Pointer to current fe_space and mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p(fe_space_p->Mesh());
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // II : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
#if SOLUTION
  // Coefficient functions used in class ReactionDiffusionElementMatrixProvider
  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0 + x.dot(x); });
  auto gamma = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  // Initialize element matrix builder
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                                      decltype(gamma)>
      elmat_builder(fe_space_p, alpha, gamma);
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return A;
};

/* SAM_LISTING_END_2 */

}  // namespace BoundaryWave
