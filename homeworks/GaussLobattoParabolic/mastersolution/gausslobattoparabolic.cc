/**
 * @file gausslobattoparabolic.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "gausslobattoparabolic.h"

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <functional>
#include <memory>
#include <utility>

namespace GaussLobattoParabolic {

/* SAM_LISTING_BEGIN_1 */
lf::assemble::COOMatrix<double> initMbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space) {
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  // Diffusion coefficient =0, reaction coefficient = 1
  lf::mesh::utils::MeshFunctionConstant alpha(0.0), gamma(1.0);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider entity_matrix_provider(
      fe_space, alpha, gamma);
  // Compute mass matrix for full finite element space
  lf::assemble::COOMatrix<double> M =
      lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>,
                                          decltype(entity_matrix_provider)>(
          0, dofh, entity_matrix_provider);
  // Find mesh nodes on the boundary
  const lf::mesh::utils::CodimMeshDataSet<bool> bd_flags =
      lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 2);
  // Predicate for selecting matrix rows induced by test functions associated
  // with nodes on the boundary
  auto pred = [&bd_flags, &dofh](int i, int j) {
    return bd_flags(dofh.Entity(i));
  };
  // Set the corresponding triplets to zero using LehrFEM++ helper function
  M.setZero(pred);

  return M;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
lf::assemble::COOMatrix<double> initAbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space) {
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  // Diffusion coefficient =1, reaction coefficient = 0
  lf::mesh::utils::MeshFunctionConstant alpha(1.0), gamma(0.0);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider entity_matrix_provider(
      fe_space, alpha, gamma);
  // Compute Galerkin matrix for -Laplacian on full FE space
  lf::assemble::COOMatrix<double> A =
      lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>,
                                          decltype(entity_matrix_provider)>(
          0, dofh, entity_matrix_provider);
  // Find mesh nodes on the boundary
  const lf::mesh::utils::CodimMeshDataSet<bool> bd_flags =
      lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 2);
  // Predicate for selecting matrix rows induced by test functions associated
  // with nodes on the boundary
  auto pred = [&bd_flags, &dofh](int i, int j) {
    return bd_flags(dofh.Entity(i));
  };
  // Set the corresponding triplets to zero using LehrFEM++ helper function
  A.setZero(pred);
  // Set "boundary block" to the identity matrix
  for (int i = 0; i < dofh.NumDofs(); ++i) {
    if (bd_flags(dofh.Entity(i))) A.AddToEntry(i, i, 1.0);
  }

  return A;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
RHSProvider::RHSProvider(const lf::assemble::DofHandler& dofh,
                         std::function<double(double)> g)
    : g_(std::move(g)) {
  // Finde nodes on the boundary
  const lf::mesh::utils::CodimMeshDataSet<bool> bd_flags =
      lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
  int N = dofh.NumDofs();
  // Initialize the fixed vector, components for degrees of freedom associated
  // with nodes on the boundary are set to 1, all other to 0
  zero_one_ = Eigen::VectorXd(N);
  for (int i = 0; i < N; ++i) {
    zero_one_(i) = bd_flags(dofh.Entity(i)) ? 1.0 : 0.0;
  }
}

Eigen::VectorXd RHSProvider::operator()(double t) const {
  // Just rescale the stored vector
  return g_(t) * zero_one_;
}
/* SAM_LISTING_END_3 */

}  // namespace GaussLobattoParabolic
