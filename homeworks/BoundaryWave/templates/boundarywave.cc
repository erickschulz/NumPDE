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
  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
  return A;
};

/* SAM_LISTING_END_2 */

}  // namespace BoundaryWave
