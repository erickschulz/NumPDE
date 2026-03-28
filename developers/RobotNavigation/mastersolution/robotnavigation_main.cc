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
  // Load mesh
  auto mesh_factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/room.msh");
  auto mesh_p = reader.mesh();
  std::cout << "Mesh: " << mesh_p->NumEntities(0) << " cells, "
            << mesh_p->NumEntities(2) << " nodes\n";

  // Solve Poisson BVP
  Eigen::VectorXd sol = RobotNavigation::solvePoissonBVP(reader);
  std::cout << "Solution range: [" << sol.minCoeff() << ", "
            << sol.maxCoeff() << "]\n";

  // Set up FE spaces
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(
          mesh_p);
  // Vector DOF handler: 2 DOFs per node (gradient components)
  lf::assemble::UniformFEDofHandler vec_dofh(
      mesh_p, {{lf::base::RefEl::kPoint(), 2},
               {lf::base::RefEl::kSegment(), 0},
               {lf::base::RefEl::kTria(), 0},
               {lf::base::RefEl::kQuad(), 0}});

  // Solve gradient projection
  Eigen::VectorXd grad =
      RobotNavigation::solveGradientProjection(fe_space, sol,
                                               vec_dofh);
  std::cout << "Gradient projection solved\n";

  // Write nodal solution and navigation field (-grad u)
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  {
    std::ofstream ofs("solution_data.csv");
    ofs << "x,y,u,gradb_x,gradb_y\n";
    for (const lf::mesh::Entity *node : mesh_p->Entities(2)) {
      auto scal_idx = dofh.GlobalDofIndices(*node)[0];
      auto vec_idx = vec_dofh.GlobalDofIndices(*node);
      auto coords = lf::geometry::Corners(*(node->Geometry()));
      ofs << coords(0, 0) << "," << coords(1, 0) << ","
          << sol(scal_idx) << ","
          << -grad(vec_idx[0]) << "," << -grad(vec_idx[1])
          << "\n";
    }
  }
  std::cout << "Wrote solution_data.csv\n";

  // Trace robot path from a starting point
  Eigen::Vector2d start(3.75, 3.0);  // between shelves 1 and 2
  Eigen::MatrixXd path = RobotNavigation::tracePath(
      mesh_p, dofh, sol, vec_dofh, grad,
      start, 0.01, 10000, 0.1);
  std::cout << "Path: " << path.rows() << " steps\n";

  {
    std::ofstream ofs("path_data.csv");
    ofs << "x,y\n";
    for (Eigen::Index i = 0; i < path.rows(); ++i) {
      ofs << path(i, 0) << "," << path(i, 1) << "\n";
    }
  }
  std::cout << "Wrote path_data.csv\n";
  return 0;
}
