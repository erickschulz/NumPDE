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

    //====================
    // Your code goes here
    // Compute the number N of DOFs and the L2-error,
    // and replace the lines below:
    int N = 1;
    double l2_error = 1.0;
    // If you want, you can use the function write_vtk(...),
    // defined in this file, to plot your solution.
    //====================

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

  //====================
  // Your code goes here
  // Replace the two variables below:
  int threshold = 0.0;  // Threshold computed by findCFLthreshold(...)
  int cfl_thres = 0.0;  // Threshold obtained from CLF using computeHmin(mesh)
                        //====================

  std::cout << "Threshold for level " << level << " is " << threshold
            << " | Threshold from CFL is " << cfl_thres << std::endl;
  /* SAM_LISTING_END_2 */

  return 0;
}
