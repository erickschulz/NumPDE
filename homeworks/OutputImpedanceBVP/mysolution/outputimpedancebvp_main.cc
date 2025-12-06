/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
#include <string>

#include "evalclass.h"
#include "outputimpedancebvp.h"

using namespace OutputImpedanceBVP;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "*** OutputImpedanceBVP ****" << std::endl;

  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/OutputImpedanceBVP.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Dirichlet boundary conditions
  Eigen::Vector2d g;
  g << 1.0, 3.0;

  // Solve BVP
  Eigen::VectorXd discrete_solution = solveImpedanceBVP(fe_space_p, g);

  // Evaluate functional
  Eigen::Vector2d d;
  d << 1.0, 2.0;
  double func_val =
      computeBoundaryOutputFunctional(discrete_solution, fe_space_p, d);
  std::cout << "Value of boundary functional : " << func_val << std::endl;

  // VERIFICATONS
  // Basis vectors for 2D Euclidean space
  Eigen::Vector2d e0, e1;
  e0 << 1.0, 0.0;
  e1 << 0.0, 1.0;
  Eigen::VectorXd approx_sol_e0 = solveImpedanceBVP(fe_space_p, e0);
  Eigen::VectorXd approx_sol_e1 = solveImpedanceBVP(fe_space_p, e1);
  double error =
      (discrete_solution - g(0) * approx_sol_e0 - g(1) * approx_sol_e1).norm();
  std::cout << "Linearity mismatch = " << error << std::endl;

  // Testing evalclass operator
  EvalResponse F_evaluator(fe_space_p);
  std::cout
      << "Value of boundary functional as computed by EvalResponse class : "
      << F_evaluator(g, d) << std::endl;

  // Output results to csv file
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::string errors_file_name = "discrete_solution.csv";
  std::ofstream file(errors_file_name.c_str());
  if (file.is_open()) {
    file << discrete_solution.format(CSVFormat);
  }

  // Output results to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "OutputImpedanceBVP_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_solution[global_idx];
  };
  vtk_writer.WritePointData("OutputImpedanceBVP_solution", *nodal_data);
  /* SAM_LISTING_END_1 */
  std::cout << "\n The OutputImpedanceBVP_solution was written to:"
            << std::endl;
  std::cout << ">> OutputImpedanceBVP_solution.vtk.vtk\n" << std::endl;
}
