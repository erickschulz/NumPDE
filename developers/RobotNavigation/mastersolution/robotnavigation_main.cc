/**
 * @file robotnavigation_main.cc
 * @brief NPDE homework RobotNavigation code
 * @author Erick Schulz
 * @date 28.03.2026
 * @copyright Developed at ETH Zurich
 */

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <fstream>
#include <iostream>

#include "robotnavigation.h"

int main() {
  // Step 1: Load mesh
  std::cout << "Loading mesh...\n";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/room.msh");
  auto mesh_p = reader.mesh();
  std::cout << "  " << mesh_p->NumEntities(0) << " cells, "
            << mesh_p->NumEntities(2) << " nodes\n";

  // Step 2: Solve Poisson BVP
  std::cout << "Solving Poisson BVP...\n";
  Eigen::VectorXd sol = RobotNavigation::solvePoissonBVP(reader);

  // Step 3: Solve gradient projection
  std::cout << "Solving gradient projection...\n";
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  lf::assemble::UniformFEDofHandler vec_dofh(mesh_p,
                                             {{lf::base::RefEl::kPoint(), 2},
                                              {lf::base::RefEl::kSegment(), 0},
                                              {lf::base::RefEl::kTria(), 0},
                                              {lf::base::RefEl::kQuad(), 0}});
  Eigen::VectorXd grad =
      RobotNavigation::solveGradientProjection(fe_space, sol, vec_dofh);

  // Step 4: Write solution data
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  std::ofstream sol_file("solution_data.csv");
  sol_file << "x,y,u,gradb_x,gradb_y\n";
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    auto scal_idx = dofh.GlobalDofIndices(*node)[0];
    auto vec_idx = vec_dofh.GlobalDofIndices(*node);
    auto coords = lf::geometry::Corners(*(node->Geometry()));
    sol_file << coords(0, 0) << "," << coords(1, 0) << "," << sol(scal_idx)
             << "," << -grad(vec_idx[0]) << "," << -grad(vec_idx[1]) << "\n";
  }
  std::cout << "  Wrote solution_data.csv\n";

  // Step 5: Trace robot path
  std::cout << "Tracing robot path...\n";
  Eigen::Vector2d start(3.75, 3.0);
  Eigen::MatrixXd path = RobotNavigation::integrateRobotPath(
      mesh_p, dofh, sol, vec_dofh, grad, start, 0.01, 10000, 0.1);
  std::ofstream path_file("path_data.csv");
  path_file << "x,y\n";
  for (Eigen::Index i = 0; i < path.rows(); ++i) {
    path_file << path(i, 0) << "," << path(i, 1) << "\n";
  }
  std::cout << "  Wrote path_data.csv (" << path.rows() << " steps)\n";
  return 0;
}
