/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"

// Eigen includes
#include <Eigen/Dense>
#include <array>

namespace LFPPDofHandling {

/**
 * @brief Compute the internal number of dofs per entity
 * @param dofhandler Dofhandler for the FE space
 * @return array with dofs per entity, index = codim
 *
 * The number of dofs per entity doesn't include the number of
 * dofs of the subordinate entities. E.g. will cells and edges in a S_1^0
 * FE space have 0 dofs.
 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler& dofhandler);

/**
 * @brief Count number of dofs on the boundary
 * @param dofhandler Dofhandler for the FE space
 * @return Number of dofs on the boundary
 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler& dofhandler);

/**
 * @brief Integrate a function over a linear FE space
 * @param dofhandler Dofhandler for the linear FE space
 * @param mu Coefficients of the basis functions, indices according to global
 *           indices of the dofs
 * @return Approximated integral
 */
double integrateLinearFEFunction(const lf::assemble::DofHandler& dofhandler,
                                 const Eigen::VectorXd& mu);

/**
 * @brief Integrate a function over a quadratic FE space
 * @param dofhandler Dofhandler for the quadratic FE space
 * @param mu Coefficients of the basis functions, indices according to global
 *           indices of the dofs
 * @return Approximated integral
 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler& dofhandler,
                                    const Eigen::VectorXd& mu);

/**
 * @brief Convert coefficient vector of linear FE space to quadratic FE space
 * @param dofh_Linear_FE Dofhandler to the linear FE space
 * @param dofh_Quadratic_FE Dofhandler to the quadratic FE space
 * @param Coefficient vector of the linear FE space
 * @return Coefficient vector of the quadratic FE space
 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu);
Eigen::VectorXd convertDOFsLinearQuadratic_alt(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu);

/**
 * @brief Evaluate the function f in the dofs and put them at the index given by
 * the dof
 * @param f Operator allowing ()-evaluation, taking an Eigen::Vector2d to a
 * double
 * @param dofh Dofhandler of indicated type
 * @return Coefficent vector of size dofh.NumDofs()
 * @note Only dofs on **nodes** and **edges** are taken into account!
 *       Not suitable for spaces S_p^0 with p > 2!
 */
template <typename F>
Eigen::VectorXd buildCoefVector(const F& f,
                                const lf::assemble::UniformFEDofHandler& dofh) {
  Eigen::VectorXd mu(dofh.NumDofs());
#define computing_mu_version 2
#if computing_mu_version == 1
  // Version 1: Iterate over nodes and edges, get corresponding dof then
  for (const auto& node : mesh->Entities(2)) {
    // coordinates of nodes
    const Eigen::Vector2d coords = lf::geometry::Corners(*node.Geometry());
    // global index
    const unsigned node_idx = dofh.InteriorGlobalDofIndices(node)[0];
    mu(node_idx) = f(coords);
  }
  for (const auto& edge : mesh->Entities(1)) {
    // coordinates of nodes
    const Eigen::Matrix2d endpoints = lf::geometry::Corners(*edge.Geometry());
    const Eigen::Vector2d midpoint =
        0.5 * (endpoints.col(0) + endpoints.col(1));
    // global index
    const unsigned edge_idx = dofh.InteriorGlobalDofIndices(edge)[0];
    mu(edge_idx) = f(midpoint);
  }
#else
  // Version 2: Iterate over dofs and get corresponding entity
  for (std::size_t dofnum = 0; dofnum < dofh.NumDofs(); ++dofnum) {
    // get associated entity
    if (dofh.Entity(dofnum).RefEl() == lf::base::RefEl::kPoint()) {
      // we're dealing with a point
      const Eigen::Vector2d coords =
          lf::geometry::Corners(*dofh.Entity(dofnum).Geometry());
      mu(dofnum) = f(coords);
    } else if (dofh.Entity(dofnum).RefEl() == lf::base::RefEl::kSegment()) {
      // we're dealing with an edge, the dof is located in its middle
      const Eigen::Matrix2d endpoints =
          lf::geometry::Corners(*dofh.Entity(dofnum).Geometry());
      const Eigen::Vector2d midpoint =
          0.5 * (endpoints.col(0) + endpoints.col(1));
      mu(dofnum) = f(midpoint);
    }
  }
#endif

  return mu;
}

}  // namespace LFPPDofHandling
