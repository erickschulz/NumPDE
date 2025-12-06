/**
 This homework problem consists of reading a simple, gmesh generated, mesh on
 the unit square and solving a simple reaction diffusion system using LehrFEM++
 */

#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>

#include "linfereactdiff.h"

int main() {
  const lf::base::size_type num_levels = 5;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      LinFeReactDiff::generateMeshHierarchy(num_levels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // get pointer to finest mesh used as ground truth
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      multi_mesh.getMesh(num_levels - 1);
  Eigen::VectorXd finest_sol = LinFeReactDiff::solveFE(mesh_p);
  double ground_truth_energy =
      LinFeReactDiff::computeEnergy(mesh_p, finest_sol);

  // compute error for the other meshes
  for (int i = 0; i < num_levels - 1; i++) {
    mesh_p = multi_mesh.getMesh(i);
    Eigen::VectorXd sol = LinFeReactDiff::solveFE(mesh_p);
    double energy = LinFeReactDiff::computeEnergy(mesh_p, sol);
    std::cout << "Mesh " << i + 1
              << " error: " << std::abs(energy - ground_truth_energy) << "\n";
  }
}
