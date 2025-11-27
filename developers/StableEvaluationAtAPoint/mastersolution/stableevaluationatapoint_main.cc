/**
 * @file stableevaluationatapoint_main.cc
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher, Erick Schulz & Philippe Peter
 * @date 29.11.2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "stableevaluationatapoint.h"
#include "systemcall.h"

int main(int /*argc*/, const char** /*argv*/) {
  // exact solution
  auto uExact = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };
  // Fixed evaluation point
  Eigen::Vector2d x(0.3, 0.4);
  std::cout << "Exact evaluation at (0.3,0.4) : " << uExact(x) << std::endl;

  // Number of meshes used in the error analysis:
  int N_meshes = 6;

  // Error analysis vectors:
  Eigen::VectorXd mesh_sizes(N_meshes);
  mesh_sizes.setZero();
  Eigen::VectorXd dofs(N_meshes);
  dofs.setZero();

  // Error vector used for the error analysis in exercise b)
  Eigen::VectorXd errors_potential(N_meshes);
  errors_potential.setZero();

  // Error vectors used for the error analysis in exercise h
  Eigen::VectorXd errors_direct(N_meshes);
  errors_direct.setZero();
  Eigen::VectorXd errors_stable(N_meshes);
  errors_stable.setZero();

  // iterate over meshes:
  for (int k = 0; k < N_meshes; k++) {
    // read mesh::
    std::string idx = std::to_string(k + 1);
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory),
                              "meshes/square" + idx + ".msh");
    auto mesh_p = reader.mesh();

    // Initialize fe-space and dofh
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
    dofs(k) = dofh.NumDofs();

    // Printing mesh statistics
    mesh_sizes(k) = StableEvaluationAtAPoint::MeshSize(mesh_p);
    std::cout << "square" + idx + ".msh: " << "N_dofs = " << dofs(k)
              << ", h=" << mesh_sizes(k) << std::endl;

    // Error anlysis part b) (Potentials)
    errors_potential(k) = StableEvaluationAtAPoint::PointEval(mesh_p);

    // error analysis part g/h: Compare direct vs stable point evaluation:
    auto [direct_eval, stable_eval] =
        StableEvaluationAtAPoint::ComparePointEval(fe_space, uExact, x);
    errors_direct(k) = std::abs(uExact(x) - direct_eval);
    errors_stable(k) = std::abs(uExact(x) - stable_eval);
  }

  // Compute rates of convergence:
  Eigen::VectorXd rates_potential(N_meshes - 1);
  Eigen::VectorXd rates_direct(N_meshes - 1);
  Eigen::VectorXd rates_stable(N_meshes - 1);

#if SOLUTION
  for (int k = 0; k < N_meshes - 1; ++k) {
    double log_denum = std::log(mesh_sizes(k) / mesh_sizes(k + 1));
    rates_potential(k) =
        std::log(errors_potential(k) / errors_potential(k + 1)) / log_denum;
    rates_direct(k) =
        std::log(errors_direct(k) / errors_direct(k + 1)) / log_denum;
    rates_stable(k) =
        std::log(errors_stable(k) / errors_stable(k + 1)) / log_denum;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  // Report computed errors and rates:
  std::cout << "Subtask b) Evaluation based on Potentials \n";
  std::cout << "Errors: \n" << errors_potential << "\n";
  std::cout << "Rates: \n" << rates_potential << "\n";
  std::cout << "Subtask h) Comparison of direct and stable evaluation: \n";
  std::cout << "Errors direct: \n" << errors_direct << "\n";
  std::cout << "Rates direct: \n" << rates_direct << "\n";
  std::cout << "Errors stable: \n" << errors_stable << "\n";
  std::cout << "Rates stable: \n" << rates_stable << "\n";

  // Output
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  Eigen::MatrixXd convergence_potential(N_meshes, 2);
  convergence_potential << mesh_sizes, errors_potential;

  Eigen::MatrixXd convergence_stable(N_meshes, 3);
  convergence_stable << mesh_sizes, errors_direct, errors_stable;

  std::ofstream file;
  file.open("convergence_potential.csv");
  file << "h, Error u(x) (Potential) \n";
  file << convergence_potential.format(CSVFormat);
  file.close();
  std::cout << "Generated convergence_potential.csv" << std::endl;

  file.open("convergence_stable.csv");
  file << "h, Error u(x) (Direct), Error u(x) (Stable) \n";
  file << convergence_stable.format(CSVFormat);
  file.close();
  std::cout << "Generated convergence_stable.csv" << std::endl;

  // Plot
  systemcall::execute("python3 scripts/plot_convergence_potential.py .");
  systemcall::execute("python3 scripts/plot_convergence_stable.py .");

  return 0;
}
