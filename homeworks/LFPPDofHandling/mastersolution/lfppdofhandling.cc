/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <lf/assemble/assembly_types.h>
#include <lf/base/lf_assert.h>
#include <lf/mesh/entity.h>

#include <Eigen/Dense>
#include <array>
#include <memory>
#include <stdexcept>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler& dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  // Idea: iterate over entities in the mesh and get interior number of dofs for
  // each
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (std::size_t codim = 0; codim <= 2; ++codim) {
    entityDofs[codim] = 0;
    for (const auto* el : mesh->Entities(codim)) {
      if (el->RefEl() == lf::base::RefEl::kQuad()) {
        throw "Only triangular meshes are allowed!";
      }
      entityDofs[codim] += dofhandler.NumInteriorDofs(*el);
    }
  }
  return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler& dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd\_flags(entity) == true, if the entity is on the
  // boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  // Edges and nodes can be on the boundary
  for (const auto* edge : mesh->Entities(1)) {
    if (bd_flags(*edge)) {
      no_dofs_on_bd += dofhandler.NumInteriorDofs(*edge);
    }
  }
  for (const auto* node : mesh->Entities(2)) {
    if (bd_flags(*node)) {
      no_dofs_on_bd += dofhandler.NumInteriorDofs(*node);
    }
  }
  return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
  double I = 0;
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (const auto *cell : mesh->Entities(0)) {
    if (cell->RefEl() == lf::base::RefEl::kQuad()) {
        throw "Only triangular meshes are allowed!";
    }    
    // check if we the FE space is really $\Cs_1^0$
    if (dofhandler.NumLocalDofs(*cell) != 3) {
      throw "Not a S_1^0 FE space on a triangular mesh!";
    }
    // Retrieve area of the triangle by calling
    // a LehrFEM++ utility function
    const double area = lf::geometry::Volume(*(cell->Geometry()));
    // iterate over dofs
    auto int_dofs = dofhandler.GlobalDofIndices(*cell);
    for (auto dof_idx_p = int_dofs.begin(); dof_idx_p < int_dofs.end();
         ++dof_idx_p) {
      // local integral of the basis function associated with this dof:
      // in linear Lagrangian FE, the integral over the basis functions over
      // a triangle K is: 1/3*vol(K)
      const double I_bary = 1.0 / 3.0 * area;
      // multiply by the value at the dof to get local contribution
      I += I_bary * mu(*dof_idx_p);
    }
  }
  return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler& dofhandler,
                                    const Eigen::VectorXd& mu) {
  double I = 0;
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (const auto* cell : mesh->Entities(0)) {
    // check if we the FE space is really $\Cs_2^0$
    if (dofhandler.NumLocalDofs(*cell) != 6) {
      throw "Not a S_2^0 FE space!";
    }
    const double weight = 1.0 / 3.0 * lf::geometry::Volume(*(cell->Geometry()));
    // iterate over dofs
    auto int_dofs = dofhandler.GlobalDofIndices(*cell);
    // The integrated basis functions associated with the nodes are 0:
    //  $\int_K b_K^j dx = 0$ for $j = 1,2,3$. Skip!
    for (int l = 3; l < 6; ++l) {
      // The integrated basis functions associated with the edges are:
      // $\int_K b_K^j dx = |K|/3$ for $j = 4,5,6$
      // multiply by the value at the dof to get local contribution
      I += (weight * mu(int_dofs[l]));
    }
  }
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();  // get the mesh
  // coefficient vector for returning the result
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());
  // Play safe: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();
  for (const auto* cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    if (dofh_Linear_FE.NumLocalDofs(*cell) != 3 ||
        dofh_Quadratic_FE.NumLocalDofs(*cell) != 6) {
      throw std::runtime_error(
          "dofh_Linear_FE must have 3 dofs per cell and dofh_Quadratic_FE 6!");
    }
    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin\_dofs will have size 3 for the 3 dofs on the nodes and
    // quad\_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    std::span<const lf::assemble::gdof_idx_t> lin_dofs =
        dofh_Linear_FE.GlobalDofIndices(*cell);
    std::span<const lf::assemble::gdof_idx_t> quad_dofs =
        dofh_Quadratic_FE.GlobalDofIndices(*cell);
    for (std::size_t l = 0; l <= 2; ++l) {
      // Let $p_1,p_2,p_3$ be the nodes of the triangle and $p_4,p_5,p_6$ the
      // edge midpoints. Let $\lambda_1,\lambda_2,\lambda_3$ be the local basis
      // functions of $\Cs_1^0$ and $b_1,..,b_6$ for $\Cs_2^0$. The values of
      // $\lambda_i$ in the interpolation nodes $p_1,..,p_6$ are the
      // coefficients of writing $\lambda_i$ as lin. comb of $b_1,..b_6$.
      // From 2-9.a we know $\lambda_l(p_l) = 1$; $l = 1,2,3$.
      // Hence we simply \textbf{copy} the coefficient of $\lambda_l$ to $b_l$
      // for $l=1,2,3$.
      zeta(quad_dofs[l]) = mu(lin_dofs[l]);
      // And $\lambda_l(p_l+3) = 0.5$, $\lambda_{(l+1) mod 3}(p_l+3) = 0.5$;
      // $l =1,2,3$. Hence we copy 0.5 of the respective coefficients!
      zeta(quad_dofs[l + 3]) =
          0.5 * (mu(lin_dofs[l]) + mu(lin_dofs[(l + 1) % 3]));
    }
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_7 */
Eigen::VectorXd convertDOFsLinearQuadratic_alt(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu) {
  LF_ASSERT_MSG(dofh_Linear_FE.Mesh() == dofh_Quadratic_FE.Mesh(),
                "Underlying meshes must be the same for both DOF handlers!");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh_Linear_FE.Mesh();
  // coefficient vector for returning the result
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());
  // Visit all nodes of the mesh and copy dof values
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    // Obtain global index numbers of the GSFs at the node for both finite
    // element spaces
    std::span<const lf::assemble::gdof_idx_t> lin_dofs =
        dofh_Linear_FE.InteriorGlobalDofIndices(*node);
    std::span<const lf::assemble::gdof_idx_t> quad_dofs =
        dofh_Quadratic_FE.InteriorGlobalDofIndices(*node);
    LF_ASSERT_MSG(lin_dofs.size() == 1, "One dof per node expected");
    LF_ASSERT_MSG(quad_dofs.size() == 1, "One dof per node expected");
    // Just copy the coefficient
    zeta[quad_dofs[0]] = mu[lin_dofs[0]];
  }
  // Run through all edges of the mesh
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    // Obtain pointers to endpoints of edge
    std::span<const lf::mesh::Entity* const> endpoints{edge->SubEntities(1)};
    LF_ASSERT_MSG(endpoints.size() == 2, "Edge must have two endpoints");
    // Obtain indices of GSFs for linear FE spaces associated with endpoints
    std::span<const lf::assemble::gdof_idx_t> lindof_p0 =
        dofh_Linear_FE.InteriorGlobalDofIndices(*endpoints[0]);
    LF_ASSERT_MSG(lindof_p0.size() == 1, "Onle one dof per vertex allowed");
    std::span<const lf::assemble::gdof_idx_t> lindof_p1 =
        dofh_Linear_FE.InteriorGlobalDofIndices(*endpoints[1]);
    LF_ASSERT_MSG(lindof_p1.size() == 1, "Onle one dof per vertex allowed");
    std::span<const lf::assemble::gdof_idx_t> quad_dofs =
        dofh_Quadratic_FE.InteriorGlobalDofIndices(*edge);
    LF_ASSERT_MSG(quad_dofs.size() == 1,
                  "Only one dof associated with an edge");
    // Set value of edge dof to average of adjacent vertex dofs
    zeta[quad_dofs[0]] = 0.5 * (mu[lindof_p0[0]] + mu[lindof_p1[0]]);
  }
  return zeta;
}
/* SAM_LISTING_END_7 */

}  // namespace LFPPDofHandling
