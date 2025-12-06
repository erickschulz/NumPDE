/**
 * @file radauthreetimestepping_main.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

#include "radauthreetimestepping.h"
#include "radauthreetimesteppingode.h"

using namespace RadauThreeTimestepping;

int main(int /*argc*/, char** /*argv*/) {
  /* Solving the ODE problem */
  // This function prints to the terminal the convergence rates and average rate
  // of a convergence study performed for the ODE (d/dt)y = -y.
  testConvergenceTwoStageRadauLinScalODE();

  /* Solving the parabolic heat equation */
  // Create a Lehrfem++ square tensor product mesh
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(50)
      .setNumYCells(50);
  auto mesh_p = builder.Build();

  /* SAM_LISTING_BEGIN_1 */
  // Generate the linear lagrange FE data
  // Finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Solve heat evolution with zero initial and boundary conditions
  double final_time = 1.0;
  unsigned int m = 50;
  Eigen::VectorXd discrete_heat_solution =
      solveHeatEvolution(dofh, m, final_time);
  LF_ASSERT_MSG(
      discrete_heat_solution.size() == N_dofs,
      "Size of discrete solution and dimension of FE space mismatch.");

  // Output results to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "discrete_heat_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_heat_solution[global_idx];
  };
  vtk_writer.WritePointData("discrete_heat_solution", *nodal_data);
  /* SAM_LISTING_END_1 */
  std::cout << "\n The discrete_heat_solution was written to:" << std::endl;
  std::cout << ">> discrete_heat_solution.vtk\n" << std::endl;

  return 0;
}
