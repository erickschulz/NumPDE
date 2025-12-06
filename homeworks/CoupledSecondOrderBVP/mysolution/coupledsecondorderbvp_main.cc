/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP
 * @author Erick Schulz
 * @date 13/11/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include "coupledsecondorderbvp.h"

using namespace CoupledSecondOrderBVP;

int main(int /*argc*/, const char** /*argv*/) {
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/hex1.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space = std::make_shared<FeSpaceLagrangeO2<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  /* Solve the coupled boundary value problem */
  double gamma = 1.0;  // reaction coefficientS
  // Right-hand side source function f
  auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::cos(x.norm()); });
  Eigen::VectorXd sol_vec = solveCoupledBVP(fe_space, gamma, f);

  /* Output results to vtk file */
  // We store data by keeping only the coefficients of nodal basis functions
  // In that sense, we are plotting the values of the solution at the vertices
  lf::io::VtkWriter vtk_writer(mesh_p, "CoupledSecondOrderBVP_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    if (dofh.Entity(global_idx).RefEl() == lf::base::RefElType::kPoint) {
      nodal_data->operator()(dofh.Entity(global_idx)) = sol_vec[global_idx];
    }
  };
  vtk_writer.WritePointData("CoupledSecondOrderBVP_solution", *nodal_data);
  /* SAM_LISTING_END_1 */
  std::cout << "\n The solution vector was written to:" << std::endl;
  std::cout << ">> CoupledSecondOrderBVP_solution.vtk\n" << std::endl;
}
