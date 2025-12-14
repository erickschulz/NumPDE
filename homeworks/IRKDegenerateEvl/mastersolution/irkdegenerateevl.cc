/**
 * @file irkdegenerateevl.cc
 * @brief NPDE homework IRKDegenerateEvl code
 * @author R. Hiptmair
 * @date June 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#include "irkdegenerateevl.h"

#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/fe/fe_tools.h>
#include <lf/fe/loc_comp_ellbvp.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/special_entity_sets.h>

#include <cstddef>
#include <memory>
#include <vector>

namespace IRKDegenerateEvl {
/* SAM_LISTING_BEGIN_1 */
lf ::assemble::COOMatrix<double> buildM(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fes_p) {
  // Extract mesh
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p{fes_p->Mesh()};
  // Fetch DofHandler
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  const size_t N = dofh.NumDofs();
  // Sparse matrix in triplet format
  lf::assemble::COOMatrix<double> M_coo(N, N);
  // Flag edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Set up ENTITY\_MATRIX\_PROVIDER
  lf::mesh::utils::MeshFunctionGlobal mf_one{
      [](Eigen::Vector2d /*x*/) -> double { return 1.0; }};
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(mf_one),
                                      decltype(bd_flags)>
      memp(fes_p, mf_one, bd_flags);
  // Assembly of matrix M
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, memp, M_coo);
  return M_coo;
}
/* SAM_LISTING_END_1 */

lf::assemble::COOMatrix<double> buildA(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fes_p) {
  // Extract mesh
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p{fes_p->Mesh()};
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  // Fetch DofHandler
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  const size_t N = dofh.NumDofs();
  // Set up ENTITY_MATRIX_PROVIDER
  lf::mesh::utils::MeshFunctionGlobal mf_one{
      [](Eigen::Vector2d /*x*/) -> double { return 1.0; }};
  lf::fe::DiffusionElementMatrixProvider<double, decltype(mf_one)> demp(fes_p,
                                                                        mf_one);
  // Assembly of matrix A
  lf::assemble::COOMatrix<double> A_coo(N, N);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, demp, A_coo);
  return A_coo;
}

}  // namespace IRKDegenerateEvl
