/** @file
 * @brief NPDE BoundaryWave
 * @author Erick Schulz
 * @date 24/07/2019
 * @copyright Developed at ETH Zurich
 */

#include <fstream>

#include "boundarywave.h"

using namespace BoundaryWave;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "*** BoundaryWave ***" << std::endl;
  std::cout << "A mixed elliptic-hyperbolic linear evolution problem "
            << std::endl;

  /* PARAMETERS */
  double T = 1.0;
  unsigned int N = 100;

  /* TOOLS AND DATA */
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/BoundaryWave.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>
  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};

  /* GENERATE INITIAL CONDITIONS */
  auto u0 = [](const Eigen::Vector2d& x) -> double {
    return x[0] + x[1] * x[1];
  };
  auto v0 = [](const Eigen::Vector2d& x) -> double {
    return 3.0 * x[0] + x[1];
  };

  /* SOLVE BOUNDARY VALUE PROBLEM */
  Eigen::VectorXd discrete_solution =
      solveBoundaryWave(fe_space_p, u0, v0, T, N);

  /* OUTPUT RESULTS */
  // Output results to csv file
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::string errors_file_name = "BoundaryWave_solution.csv";
  std::ofstream file(errors_file_name.c_str());
  if (file.is_open()) {
    file << discrete_solution.format(CSVFormat);
  }

  // Output results to vtk file
  std::string vtk_file_name = "BoundaryWave_solution.vtk";
  lf::io::VtkWriter vtk_writer(mesh_p, vtk_file_name);
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < discrete_solution.rows();
       global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_solution[global_idx];
  };
  vtk_writer.WritePointData("BoundaryWave_solution", *nodal_data);

  std::cout << "\n The BoundaryWave_solution was written to:" << std::endl;
  std::cout << ">> " << vtk_file_name << std::endl;
}
