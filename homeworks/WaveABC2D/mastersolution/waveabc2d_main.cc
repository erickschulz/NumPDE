/** @file
 * @brief NPDE WaveABC2D
 * @author Erick Schulz
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
#include <memory>
#include <string>

#include "waveabc2d.h"
// Eigen includes
#include <Eigen/Core>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

using namespace WaveABC2D;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "\n" << std::endl;
  std::cout << "PROBLEM - WaveABC2D" << std::endl;

  testConvergenceScalarImplicitTimestepping();

  // Load mesh into a Lehrfem++ object
  std::cout << "Loading mesh..." << std::endl;
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/unitsquare2.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  auto mu0 = [](const Eigen::Vector2d& x) -> double {
    return std::sin(x.norm());
  };
  auto nu0 = [](const Eigen::Vector2d& x) -> double { return std::cos(x(1)); };
  auto rho = [](Eigen::Vector2d) -> double { return 1.0; };

  WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)> stepper(
      fe_space_p, rho, 250, 1.0);
  Eigen::VectorXd discrete_solution = stepper.solveWaveABC2D(mu0, nu0);

  double discrete_energy = stepper.energies();

  // Output results to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "WaveABC2D_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_solution[global_idx];
  };
  vtk_writer.WritePointData("WaveABC2D_solution", *nodal_data);

  std::cout << "\nThe WaveABC2D_solution was written to:" << std::endl;
  std::cout << "WaveABC2D_solution.vtk\n" << std::endl;

  std::cout << "The discrete energies E^(k) : " << discrete_energy << std::endl;
  return 0;
}  // main
