/**
 * @file
 * @brief  NPDE UnstableBVP. Solution of source-free heat equation and
 * computation of H1 seminorms on different triangular meshes and refinement
 * levels
 * @author Julien Gacon, Ralf Hiptmair, Amélie Loher
 * @date   March 2019
 * @copyright MIT License
 */
#include "unstablebvp.h"
// General includes
#include <fstream>
#include <iomanip>
#include <memory>
// Math includes
#include <cmath>
// Eigen
#include <Eigen/Core>
// Lehrfempp
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

int main() {
  // Define the number of refinement levels we want for our mesh
  const int reflevels = 7;

  // Do the computations for all three possible mesh locations:
  // Above x2 = 0 (top), intersecting it (center), or below it (bottom)
  // Store the data in a Eigen Matrix for nice display later on
  Eigen::MatrixXd h1_seminorms(reflevels + 1, 3), h1_diffs(reflevels + 1, 3);

  // 0 for top (-> first column), 1 for center, 2 for bottom
  int type_idx = 0;
  for (auto mesh_type : {"top", "center", "bottom"}) {
    // Get a hierachy of refined meshes
    std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
        UnstableBVP::createMeshHierarchy(reflevels, mesh_type);
    lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

    // Number of levels
    const int L = multi_mesh.NumLevels();

    // Computing the difference u_L - u_k, k=0...L-1 is hard, since we would
    // the functions are defined on different meshes.
    // To avoid this caveat we can approximate using the triangle inequality:
    // abs(norm(u_L) - norm(u_k)) <= norm(u_L - u_k)
    // This should also display the right convergence properties.

    // H1 seminorm of most refined solution (level L-1)
    double h1_uL =
        UnstableBVP::solveTemperatureDistribution(multi_mesh.getMesh(L - 1));

    // Store H1 seminorms and the difference of the seminorms to uL
    h1_seminorms(reflevels, type_idx) = h1_uL;
    h1_diffs(reflevels, type_idx) = 0;

    // Compute H1 seminorm for all levels 0..L-2
    for (int level = 0; level < L - 1; ++level) {
      // Get the mesh pointer
      std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);

      // Get the seminorm
      double h1 = UnstableBVP::solveTemperatureDistribution(mesh_p);

      // Store
      h1_seminorms(level, type_idx) = h1;
      h1_diffs(level, type_idx) = std::abs(h1_uL - h1);
    }
    ++type_idx;
  }

  // Print to terminal
  std::cout << std::left << std::setfill('-');
  std::cout << std::setw(39) << "#  -- Above x2 = 0 " << std::setw(30)
            << "|-- Intersecting x2 = 0 " << std::setw(30)
            << "|-- Below x2 = 0 "
            << "\n"
            << std::right << std::setfill(' ');
  std::cout << "#" << std::setw(8) << "Level" << std::setw(10) << "|u_k|"
            << std::setw(20) << "||u_L| - |u_k||" << std::setw(10) << "|u_k|"
            << std::setw(20) << "||u_L| - |u_k||" << std::setw(10) << "|u_k|"
            << std::setw(20) << "||u_L| - |u_k||"
            << "\n";

  for (int l = 0; l < reflevels; ++l) {
    std::cout << std::setw(9) << l << std::setw(10) << h1_seminorms(l, 0)
              << std::setw(20) << h1_diffs(l, 0) << std::setw(10)
              << h1_seminorms(l, 1) << std::setw(20) << h1_diffs(l, 1)
              << std::setw(10) << h1_seminorms(l, 2) << std::setw(20)
              << h1_diffs(l, 2) << "\n";
  }
  std::cout << std::setw(9) << reflevels << std::setw(10)
            << h1_seminorms(reflevels, 0) << std::setw(20) << 0 << std::setw(10)
            << h1_seminorms(reflevels, 1) << std::setw(20) << 0 << std::setw(10)
            << h1_seminorms(reflevels, 2) << std::setw(20) << 0 << "\n";

  return 0;
}
