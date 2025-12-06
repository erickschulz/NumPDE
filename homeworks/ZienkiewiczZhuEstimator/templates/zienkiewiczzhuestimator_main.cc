/** @file
 * @brief NPDE ZienkiewiczZhuEstimator
 * @author Erick Schulz
 * @date 25/07/2019
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "zienkiewiczzhuestimator.h"
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

using namespace ZienkiewiczZhuEstimator;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "\n" << std::endl;
  std::cout << "PROBLEM - ZienkiewiczZhuEstimator " << std::endl;
  progress_bar progress{std::clog, 70u, "Computing"};
  double progress_pourcentage;
  // Tools and data
  int N_meshes = 4;             // num. of meshes
  Eigen::VectorXd approx_sol;   // basis ceoff expansion of scalar approx sol
  Eigen::VectorXd approx_grad;  // basis ceoff expansion of approx grad
  Eigen::VectorXd L2errors(N_meshes);     // L2 errors of scalar approx sol
  Eigen::VectorXd H1errors(N_meshes);     // H1 errors of scalar approx sol
  Eigen::VectorXd errors_grad(N_meshes);  // deviation of grad (delta error)
  Eigen::VectorXd errors_diff(N_meshes);  // error difference (epsilon error)
  Eigen::VectorXd mesh_sizes(N_meshes);
  Eigen::VectorXd interpolated_uExact;
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p;
  std::shared_ptr<const lf::mesh::Mesh> mesh_p;

  // Analytic solution
  auto uExact = [](coord_t x) -> double {
    return (0.2 / (M_PI * M_PI)) * sin(M_PI * x[0]) * sin(2 * M_PI * x[1]);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_uExact{uExact};
  auto grad_uExact = [](coord_t x) -> Eigen::Vector2d {
    return Eigen::Vector2d(
        (0.2 / M_PI) * cos(M_PI * x[0]) * sin(2 * M_PI * x[1]),
        (0.4 / M_PI) * sin(M_PI * x[0]) * cos(2 * M_PI * x[1]));
  };
  lf::mesh::utils::MeshFunctionGlobal mf_grad_uExact{grad_uExact};

  for (int i = 0; i < N_meshes; i++) {  // for each mesh
    std::string idx_str = std::to_string(i);
    // Load mesh into a Lehrfem++ object
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory),
                                    "meshes/unitsquare" + idx_str + ".msh");
    mesh_p = reader.mesh();
    mesh_sizes[i] = getMeshSize(mesh_p);
    fe_space_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain reference to scalar dofh
    const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
    // Produce a dof handler for the vector-valued finite element space
    lf::assemble::UniformFEDofHandler vec_dofh(
        mesh_p, {{lf::base::RefEl::kPoint(), 2},
                 {lf::base::RefEl::kSegment(), 0},
                 {lf::base::RefEl::kTria(), 0},
                 {lf::base::RefEl::kQuad(), 0}});

    LF_VERIFY_MSG(2 * dofh.NumDofs() == vec_dofh.NumDofs(),
                  "Number of degrees of freedom mismatch!");

    // Solve Poisson BVP with essential BCs
    approx_sol = solveBVP(fe_space_p);
    // Compute L2 and H1 errors
    auto mf_approx_sol = lf::fe::MeshFunctionFE(fe_space_p, approx_sol);
    L2errors[i] = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::uscalfe::squaredNorm(mf_uExact - mf_approx_sol), 2));
    auto mf_approx_grad_sol =
        lf::fe::MeshFunctionGradFE(fe_space_p, approx_sol);
    H1errors[i] = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::uscalfe::squaredNorm(mf_grad_uExact - mf_approx_grad_sol),
        2));

    // Solve gradient VP and compute L2 deviation
    approx_grad = solveGradVP(fe_space_p, approx_sol, vec_dofh);
    // approx_grad = computeLumpedProjection(dofh, approx_sol, vec_dofh);
    // Compute deviation error (delta error)
    errors_grad[i] =
        computeL2Deviation(dofh, approx_sol, vec_dofh, approx_grad);

    // Compute difference error (epsilon error)
    errors_diff[i] = std::abs(H1errors[i] - errors_grad[i]);

    progress_pourcentage = ((double)i + 1.0) / N_meshes * 100.0;
    progress.write(progress_pourcentage / 100.0);
  }

  // Computing rates of convergence of approx solutions to BVP and gradient VP
  double ratesL2[N_meshes - 1];
  double ratesH1[N_meshes - 1];
  double rates_grad[N_meshes - 1];
  double rates_diff[N_meshes - 1];
  double log_denum;
  for (int k = 0; k < N_meshes - 1; k++) {
    log_denum = log(mesh_sizes[k] / mesh_sizes[k + 1]);
    ratesL2[k] = log(L2errors[k] / L2errors[k + 1]) / log_denum;
    ratesH1[k] = log(H1errors[k] / H1errors[k + 1]) / log_denum;
    rates_grad[k] = log(errors_grad[k] / errors_grad[k + 1]) / log_denum;
    rates_diff[k] = log(errors_diff[k] / errors_diff[k + 1]) / log_denum;
  }

  std::cout << "\n Done." << std::endl;
  std::cout << "Outputing results." << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "      ERRORS AND CONVERGENCE RATES FOR PROBLEM        "
            << std::endl;
  std::cout << "             ZIENKIEWICZ ZHU ESTIMATOR                   "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "    L2 errors of approximate scalar solutions to the     "
            << std::endl;
  std::cout << "   Poisson Problem with essential boundary conditions    "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| L2 error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << L2errors[k];
    if (k > 0) {
      std::cout << "\t \t|" << ratesL2[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "    H1 errors of approximate scalar solutions to the     "
            << std::endl;
  std::cout << "   Poisson Problem with essential boundary conditions    "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| H1 error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << H1errors[k];
    if (k > 0) {
      std::cout << "\t \t|" << ratesH1[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "      Convergence rates of the approximate gradients     "
            << std::endl;
  std::cout << "                     (delta error)                       "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| delta"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << errors_grad[k];
    if (k > 0) {
      std::cout << "\t\t|" << rates_grad[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "      Difference error between the H1 error of the       "
            << std::endl;
  std::cout << "   approximate solution to the BVP and the deviation     "
            << std::endl;
  std::cout << "                    (epsilon error)                      "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| epsilon"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << errors_diff[k];
    if (k > 0) {
      std::cout << "\t \t|" << rates_diff[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  // Output results for epsilon and delta errors to csv file
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::string errors_file_name_delta = "delta_errors.csv";
  std::ofstream file_delta(errors_file_name_delta.c_str());
  if (file_delta.is_open()) {
    file_delta << errors_grad.format(CSVFormat);
  }
  std::string errors_file_name_epsilon = "epsilon_errors.csv";
  std::ofstream file_epsilon(errors_file_name_epsilon.c_str());
  if (file_epsilon.is_open()) {
    file_epsilon << errors_diff.format(CSVFormat);
  }

  std::cout << "\n The delta and epsilon errors were written to:" << std::endl;
  std::cout << ">> delta_errors.csv and epsilon_errors.csv\n" << std::endl;

  // Save approximate solution in VTK format
  // Output results to vtk file
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  lf::io::VtkWriter vtk_writer(mesh_p, "ZienkiewiczZhuEstimator_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the
  // vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) = approx_sol[global_idx];
  };
  vtk_writer.WritePointData("ZienkiewiczZhuEstimator_solution", *nodal_data);

  std::cout << "\n The ZienkiewiczZhuEstimator_solution was written to:"
            << std::endl;
  std::cout << ">> ZienkiewiczZhuEstimator_solution.vtk\n" << std::endl;

}  // main
