/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include <array>
#include <iostream>
#include <memory>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lfppdofhandling.h"

int main() {
  // Create a triangular mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Create a dofhandler for the linear lagrangian FE space
  // For the S_1^0 space we have 1 dof per node, and 0 for all other entities:
  unsigned ndof_node = 1, ndof_edge = 0, ndof_tria = 0, ndof_quad = 0;
  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), ndof_node},
             {lf::base::RefEl::kSegment(), ndof_edge},
             {lf::base::RefEl::kTria(), ndof_tria},
             {lf::base::RefEl::kQuad(), ndof_quad}});

  // and for the quadratic LFES, S_2^0, edges also have 1 dof
  ndof_edge = 1;
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), ndof_node},
             {lf::base::RefEl::kSegment(), ndof_edge},
             {lf::base::RefEl::kTria(), ndof_tria},
             {lf::base::RefEl::kQuad(), ndof_quad}});

  // Subproblem c)
  std::cout << "\n -- Subproblem (c)\n";
  std::array<std::size_t, 3> entityDofs =
      LFPPDofHandling::countEntityDofs(lin_dofh);
  std::array<std::string, 3> entityNames = {"Cells", "Edges", "Nodes"};
  for (std::size_t codim = 0; codim <= 2; ++codim) {
    std::cout << entityNames[codim] << ": " << entityDofs[codim] << " dofs\n";
  }

  // Subproblem d)
  std::cout << "\n -- Subproblem (d)\n";
  std::cout << "Dofs on boundary: "
            << LFPPDofHandling::countBoundaryDofs(lin_dofh) << "\n";

  // Subproblem e)
  std::cout << "\n -- Subproblem (e)\n";
  // integrating a linear function
  auto f = [](const Eigen::Vector2d& x) { return x[0] * x[1]; };
  Eigen::VectorXd mu = LFPPDofHandling::buildCoefVector(f, lin_dofh);
  std::cout << "Integrating u..\n";
  std::cout << ".. with linear basis functions: I = "
            << LFPPDofHandling::integrateLinearFEFunction(lin_dofh, mu) << "\n";

  // Subproblem f)
  std::cout << "\n -- Subproblem (f)\n";
  // zeta is the coeffiecient for the S_2^0 space, as mu is for the S_1^0 space
  Eigen::VectorXd zeta = LFPPDofHandling::buildCoefVector(f, quad_dofh);
  std::cout << ".. with quadratic basis functions: I = "
            << LFPPDofHandling::integrateQuadraticFEFunction(quad_dofh, zeta)
            << "\n";

  // Subproblem g)
  std::cout << "\n -- Subproblem (g)\n";
  Eigen::VectorXd zeta_constructed =
      LFPPDofHandling::convertDOFsLinearQuadratic(lin_dofh, quad_dofh, mu);
  std::cout << "Relative difference of zeta constructed directly (from "
               "coordinates) and from mu: "
            << (zeta - zeta_constructed).norm() / zeta.norm() << "\n";
  return 0;
}
