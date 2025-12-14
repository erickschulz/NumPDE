#ifndef SYMPLECTIC_ASSEMBLE_HPP
#define SYMPLECTIC_ASSEMBLE_HPP
/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

// common includes
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace SymplecticTimesteppingWaves {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNC_ALPHA, typename FUNC_GAMMA, typename FUNC_BETA>
Eigen::SparseMatrix<double> assembleGalerkinMatrix(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNC_ALPHA&& alpha, FUNC_GAMMA&& gamma, FUNC_BETA&& beta) {
  Eigen::SparseMatrix<double> galMat;

  /* Creating coefficient-functions as Lehrfem++ mesh functions */
  // Coefficient-functions used in the class template
  // ReactionDiffusionElementMatrixProvider<SCALAR,DIFF_COEFF,REACTION_COEFF>
  auto alpha_mf = lf::mesh::utils::MeshFunctionGlobal(alpha);
  auto gamma_mf = lf::mesh::utils::MeshFunctionGlobal(gamma);
  // Coefficient-function used in the class template
  // MassEdgeMatrixProvider< SCALAR, COEFF, EDGESELECTOR >
  auto beta_mf = lf::mesh::utils::MeshFunctionGlobal(beta);

  /* Retrieving FE data */
  // pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Creating a predicate that will guarantee that the computations for the
  // boundary Mass matrix part of the full Galerkin matrix are carried only on
  // the edges of the mesh using the boundary flags
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge);
  };

  /* Instantiating FE Galerkin matrix */
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> galMat_COO(N_dofs, N_dofs);

  /* Initialization of local matrices builders */
  // Initialize objects taking care of local computations for volume integrals
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(alpha_mf), decltype(gamma_mf)>
      elmat_builder(fe_space_p, alpha_mf, gamma_mf);
  // Initialize objects taking care of local computations for boundary integrals
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(beta_mf),
                                      decltype(edges_predicate)>
      edgemat_builder(fe_space_p, beta_mf, edges_predicate);

  /* Assembling the Galerkin matrix */
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, galMat_COO);
  // Invoke assembly on edges by specifying co-dimension = 1
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder,
                                      galMat_COO);
  /* Convert COO matrix A into sparse matrix format using Eigen's */
  galMat = galMat_COO.makeSparse();
  return galMat;
}
/* SAM_LISTING_END_1 */

}  // namespace SymplecticTimesteppingWaves

#endif
