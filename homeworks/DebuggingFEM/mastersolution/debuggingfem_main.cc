/**
 * @file debuggingfem_main.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "locallaplaceqfe.h"
#include "qfeinterpolator.h"
#include "qfeprovidertester.h"
#include "systemcall.h"

using size_type = lf::base::size_type;

int main() {
  // read mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/square_64.msh");
  auto mesh = reader.mesh();

  // refine mesh
  const int reflevels = 4;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  lf::base::size_type L = multi_mesh.NumLevels();

  // vector holding pointers to the different element matrix providers
  std::vector<std::unique_ptr<DebuggingFEM::EntityMatrixProvider>>
      element_matrix_provider;
  element_matrix_provider.emplace_back(
      std::make_unique<DebuggingFEM::LocalLaplaceQFE1>());
  element_matrix_provider.emplace_back(
      std::make_unique<DebuggingFEM::LocalLaplaceQFE2>());
  element_matrix_provider.emplace_back(
      std::make_unique<DebuggingFEM::LocalLaplaceQFE3>());
  const int num_emp = element_matrix_provider.size();

  // vectors accumulating the computed errors
  Eigen::MatrixXd H1SMerr(L, num_emp);
  Eigen::VectorXi N(L);

  for (int level = 0; level < L; ++level) {
    // set up fespace and dof handler for the mesh at the current level
    auto mesh_p = multi_mesh.getMesh(level);
    lf::uscalfe::FeSpaceLagrangeO2<double> fespace(mesh_p);
    const auto& dofh = fespace.LocGlobMap();

    // Dimension of finite element space
    N[level] = dofh.NumDofs();

    // compute error for each element matrix provider
    for (int i = 0; i < num_emp; ++i) {
      // function to interpolateOntoQuadFE
      auto f = [](const Eigen::VectorXd& x) -> double {
        return std::exp(x(1) * x(1) + x(0) * x(0));
      };
      double energy = 0.0;
      // Matrix in triplet format holding Galerkin for LocalLaplaceQFEX matrix,
      // zero initially.
      DebuggingFEM::QFEProviderTester qfe_provider_tester(
          dofh, *element_matrix_provider[i]);
      // compute the energy and error
      energy = qfe_provider_tester.energyOfInterpolant(f);
      H1SMerr(level, i) = std::abs(23.76088 - energy);
    }
  }

  // Tabular output of the results
  std::cout << std::left << std::setw(10) << "N" << std::setw(20)
            << "Assember 1" << std::setw(20) << "Assembler 2" << std::setw(20)
            << "Assembler 3" << std::endl;
  for (int level = 0; level < L; ++level) {
    std::cout << std::left << std::setw(10) << N[level] << std::setw(20)
              << H1SMerr(level, 0) << std::setw(20) << H1SMerr(level, 1)
              << std::setw(20) << H1SMerr(level, 2) << std::endl;
  }

  // Write .csv file and plot it by a python script
  Eigen::MatrixXd data(L, num_emp + 1);
  data.col(0) = N.cast<double>();
  data.rightCols(num_emp) = H1SMerr;
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream error_file;
  error_file.open("error.csv");
  error_file << data.format(CSVFormat) << std::endl;
  error_file.close();
  std::cout << "Generated error.csv" << std::endl;
  systemcall::execute("python3 scripts/plot_error.py error.csv error.eps");

  return 0;
}
