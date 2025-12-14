/**
 * @file advectionfv2d_main.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "advectionfv2d.h"
#include "systemcall.h"

// Use this function to plot your solution
void write_vtk(const lf::assemble::DofHandler& dofh,
               const Eigen::VectorXd& solution, std::string name) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  lf::io::VtkWriter vtk_writer(mesh_p, name + ".vtk");
  auto cell_data_ref =
      lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 0);
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    int row = dofh.GlobalDofIndices(*cell)[0];
    cell_data_ref->operator()(*cell) = solution[row];
  }
  vtk_writer.WriteCellData(name, *cell_data_ref);
}

int main() {
  /* SAM_LISTING_BEGIN_1 */
  // Define velocity field beta
  // Note that the problem description requires ||B|| <= 1
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  // Functor for initial bump
  auto u0 = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d x0(0.8, 0.2);
    double d = 0.2;
    double dist = (x - x0).norm();
    if (dist < d) {
      double cos = std::cos(M_PI / (2.0 * d) * dist);
      return cos * cos;
    }
    return 0.0;
  };

  double T = 1.0;
  // Generate a mesh hierarchy
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 6)};

  std::vector<int> vector_num_cells;
  std::vector<double> vector_l2error;

  // Iterate over mesh levels starting from thrid refinement
  int num_meshes = mesh_seq_p->NumLevels();
  for (int level = 3; level < num_meshes; ++level) {
    std::cout << "Computing L2Error for level: " << level << std::endl;

    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);

    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(
        cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                   {lf::base::RefEl::kSegment(), 0},
                   {lf::base::RefEl::kTria(), 1},
                   {lf::base::RefEl::kQuad(), 1}});
    int N = cur_dofh.NumDofs();

    // Compute cell normals
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors = AdvectionFV2D::computeCellNormals(cur_dofh.Mesh());

    // Compute adjecent cells
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        std::array<const lf::mesh::Entity*, 4>>>
        adjacentCells = AdvectionFV2D::getAdjacentCellPointers(cur_dofh.Mesh());

    // Get approximate solution from simulation
    Eigen::VectorXd mu_approx = AdvectionFV2D::simulateAdvection(
        cur_dofh, beta, u0, adjacentCells, normal_vectors, T);
    write_vtk(cur_dofh, mu_approx, "approx" + std::to_string(level));

    // Get exact solution at barycenters of cells
    Eigen::VectorXd mu_exact = AdvectionFV2D::refSolution(cur_dofh, u0, T);
    write_vtk(cur_dofh, mu_exact, "exact" + std::to_string(level));

    // Compute L2 error in barycenter
    double l2_error = 0;
    for (const lf::mesh::Entity* cell : cur_mesh->Entities(0)) {
      const lf::geometry::Geometry* geo_p = cell->Geometry();
      double area = lf::geometry::Volume(*geo_p);
      int idx = cur_dofh.GlobalDofIndices(*cell)[0];
      l2_error += std::pow((mu_approx[idx] - mu_exact[idx]), 2) * area;
    }
    l2_error = std::sqrt(l2_error);

    vector_num_cells.push_back(N);
    vector_l2error.push_back(l2_error);
  }

  // Write Output file of Task 8-8.o
  std::ofstream csv_file;
  csv_file.open("advectionfv2d.csv");
  for (int i = 0; i < vector_num_cells.size(); ++i) {
    std::cout << "Cells: " << vector_num_cells.at(i)
              << " | L2Error: " << vector_l2error.at(i) << std::endl;
    csv_file << vector_num_cells.at(i) << "," << vector_l2error.at(i) << "\n";
  }
  csv_file.close();

  systemcall::execute(
      "python3 scripts/advectionfv2d.py advectionfv2d.csv "
      "solution.eps");
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_2 */
  // Task 8-8.q
  // Compute threshold for fourth refinement level
  int level = 4;
  auto mesh = mesh_seq_p->getMesh(level);

  // Create a DOF Hander for the current mesh
  const lf::assemble::UniformFEDofHandler dofh(
      mesh, {{lf::base::RefEl::kPoint(), 0},
             {lf::base::RefEl::kSegment(), 0},
             {lf::base::RefEl::kTria(), 1},
             {lf::base::RefEl::kQuad(), 1}});

  int threshold = AdvectionFV2D::findCFLthreshold(dofh, beta, T);
  int cfl_thres = int((T / AdvectionFV2D::computeHmin(mesh) + 1));

  std::cout << "Threshold for level " << level << " is " << threshold
            << " | Threshold from CFL is " << cfl_thres << std::endl;
  /* SAM_LISTING_END_2 */

  return 0;
}
