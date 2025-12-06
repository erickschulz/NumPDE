/**
 * @file
 * @brief NPDE homework ErrorEstimatesForTraces
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "systemcall.h"
#include "teelaplrobinassembly.h"

using namespace ErrorEstimatesForTraces;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "NUMPDE PROBLEM 3-5 " << std::endl;

  int N_meshes = 4;
  Eigen::MatrixXd results(N_meshes, 2);

  for (int i = 1; i <= N_meshes; i++) {  // for each mesh
    /* SAM_LISTING_BEGIN_1 */
    // Load mesh into a Lehrfem++ object
    std::string mesh_file = "meshes/hex" + std::to_string(i) + ".msh";
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
    auto mesh_p = reader.mesh();  // type shared_ptr<const lf::mesh::Mesh>

    // Finite element space
    auto fe_space = std::make_shared<linear_lagrange>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    // Dimension of finite element space
    const lf::base::size_type N_dofs(dofh.NumDofs());
    results(i - 1, 0) = N_dofs;

    // Solve the boundary value problem with Robin boundary conditions
    Eigen::VectorXd sol_vec = solveBVP(fe_space);

    // Integrate the solution sol_vec over the flagged edges
    double bd_functional_val = bdFunctionalEval(fe_space, sol_vec);

    double error = std::abs(bd_functional_val - 2.081541059732923);
    results(i - 1, 1) = error;
    /* SAM_LISTING_END_1 */

    std::cout << mesh_file;
    std::cout << "\t(Ndofs = " << N_dofs;
    std::cout << "):";
    std::cout << std::setprecision(16);
    std::cout << "\t value: " << bd_functional_val;
    std::cout << std::setprecision(16);
    std::cout << "\t error: " << error << std::endl;
  }

  // Define output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  // Write N_dofs and error to results.csv
  std::ofstream file;
  file.open("results.csv");
  file << results.format(CSVFormat);
  file.close();
  std::cout << "Generated results.csv" << std::endl;

  // Apply plot.py to results.csv
  systemcall::execute("python3 scripts/plot.py results.csv results.eps");
}
