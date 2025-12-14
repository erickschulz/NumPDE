/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
#include <string>

#include "sdirkmethodoflines.h"
#include "sdirkmethodoflines_ode.h"
#include "systemcall.h"

using namespace SDIRKMethodOfLines;

int main(int /*argc*/, char** /*argv*/) {
  /* SDIRK-2 ODE convergence */
  sdirk2ScalarODECvTest();

  std::cout << "\n************* SDIRKMethodOfLines *************" << std::endl;
  /* Solving the parabolic temperature convection cooling problem */
  // Create a Lehrfem++ square tensor product mesh
  // Obtain mesh factory
  /* std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
       std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
   // Triangular tensor product mesh
   lf::mesh::hybrid2d::TPTriagMeshBuilder builder(mesh_factory_ptr);
   // Set mesh parameters following the Builder pattern
   // Domain is the unit square
   builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
       .setTopRightCorner(Eigen::Vector2d{1, 1})
       .setNoXCells(100)
       .setNoYCells(100);
   std::shared_ptr<lf::mesh::Mesh> mesh_p{builder.Build()}; */

#if SOLUTION
  /* SAM_LISTING_BEGIN_1 */

  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/square64_bnd.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Building initial condition vector
  Eigen::VectorXd initial_temperature_vec(N_dofs);
  for (int idx = 0; idx < N_dofs; idx++) {
    // Obtain coordinates of vertex at global index idx
    auto coords = lf::geometry::Corners(*(dofh.Entity(idx).Geometry()));
    LF_ASSERT_MSG(coords.cols() == 1, "Wrong no of coords in vertex");
    initial_temperature_vec(idx) = 5.0;
  }

  // SDIRK-2 evolution of parabolic problem
  unsigned int m = 100;
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solution_pair =
      solveTemperatureEvolution(dofh, m, 1.0, initial_temperature_vec);
  Eigen::VectorXd discrete_temperature_sol = solution_pair.first;
  LF_ASSERT_MSG(
      discrete_temperature_sol.size() == N_dofs,
      "Size of discrete solution and dimension of FE space mismatch.");
  Eigen::VectorXd energies = solution_pair.second;
  LF_ASSERT_MSG(energies.size() == m + 1, "Wrong number of energie values.");
  // Define output file format for the energies
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  // Corresponding time grid for plotting
  Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(m + 1, 0.0, 1.0);

  // Write .csv file of energy vs. time
  std::ofstream file;
  file.open("energies.csv");
  file << time.transpose().format(CSVFormat) << std::endl;
  file << energies.transpose().format(CSVFormat) << std::endl;
  file.close();
  std::cout << "Generated energies.csv" << std::endl;

  // Plot from .csv file using python
  systemcall::execute(
      "python3 ms_scripts/plot_energies.py energies.csv "
      "energies.eps");
  /* SAM_LISTING_END_1 */
#else
  //====================
  // Your code goes here
  //====================
#endif

#if SOLUTION
  /* SAM_LISTING_BEGIN_2 */
  // Output results for the temperature function to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "discrete_temperature_sol.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_temperature_sol[global_idx];
  };
  vtk_writer.WritePointData("discrete_temperature_sol", *nodal_data);
  /* SAM_LISTING_END_2 */
  std::cout << "\n>>The discrete_heat_solution was written to:" << std::endl;
  std::cout << "\t discrete_temperature_sol.vtk\n" << std::endl;
#else
//====================
// Your code goes here
//====================
#endif
  return 0;
}
