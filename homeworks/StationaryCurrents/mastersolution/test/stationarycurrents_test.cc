/**
 * @ file
 * @ brief NPDE homework problem: Computation of stationary currents, Tests
 * @ author Ralf Hiptmair
 * @ date August 2020
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "../stationarycurrents.h"

#include <gtest/gtest.h>

namespace dmxbc {

TEST(StationaryCurrents, solveMixedBVP) {
  std::cout << "Test of function solveMIxedBVP" << std::endl;
  // Path to mesh .msh file
  std::string mesh_path = "meshes/bentwire0.msh";
  // Read mesh and label nodes
  auto [mesh_p, edgeids] = readMeshWithTags(mesh_path);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Distribute tags to nodes
  auto nodeids{tagNodes(mesh_p, edgeids)};
  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Coefficient function for conductivity
  auto sigma = [](Eigen::Vector2d x) -> double {
    return std::max(1.0, x[0] - 1.0);
  };
  // Solve discrete boundary value problem. Potential values at contacts are
  // passed as a list of values corresponding to "nodal ids" 0,1,... stored in
  // the nodeids array.
  Eigen::VectorXd sol_vec{solveMixedBVP(fe_space, nodeids, {0.0, 1.0}, sigma)};
  // Compute Euclidean norm of solution vector
  const double sol_norm = sol_vec.norm();
  EXPECT_NEAR(sol_norm, 4.75759, 1.0E-3)
      << "Euclidean norm of solution vector = " << sol_norm
      << ", should be  4.75759" << std::endl;
}

TEST(StationaryCurrents, stabFlux) {
  std::cout << "Test of function stabFlux" << std::endl;
  // Path to mesh .msh file
  std::string mesh_path = "meshes/bentwire0.msh";
  // Read mesh and label nodes
  auto [mesh_p, edgeids] = readMeshWithTags(mesh_path);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Distribute tags to nodes
  auto nodeids{tagNodes(mesh_p, edgeids)};
  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Coefficient function for conductivity
  auto sigma = [](Eigen::Vector2d x) -> double {
    return std::max(1.0, x[0] - 1.0);
  };
  // Solve discrete boundary value problem. Potential values at contacts are
  // passed as a list of values corresponding to "nodal ids" 0,1,... stored in
  // the nodeids array.
  Eigen::VectorXd sol_vec{solveMixedBVP(fe_space, nodeids, {0.0, 1.0}, sigma)};
  // Compute total current through one contact
  const double stab_flux = stabFlux(
      fe_space, sol_vec, sigma, [](Eigen::Vector2d x) -> Eigen::Vector2d {
        double gx = 0.0;
        if ((x[0] > 2.0) && (x[0] < 3.0)) {
          const double a = (3.0 - x[0]) * M_PI / 2;
          gx = M_PI * std::cos(a) * std::sin(a);
        }
        return Eigen::Vector2d(gx, 0.0);
      });
  EXPECT_NEAR(stab_flux, 0.190757, 1.0E-3)
      << "Total current = " << stab_flux << ", should be 0.190757" << std::endl;
}

}  // namespace dmxbc
