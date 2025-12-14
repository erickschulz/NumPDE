/**
 * @ file avgvalboundary_main.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

#include "avgvalboundary.h"

/* SAM_LISTING_BEGIN_1 */
int main() {
  // read in mesh and set up finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/square.msh");
  auto mesh = reader.mesh();
  // obtain dofh for lagrangian finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Solve test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute H1 seminorm of the solution
  double h1s_norm = AvgValBoundary::compH1seminorm(dofh, mu);
  // compute boundary functional
  auto w = [](Eigen::Vector2d x) -> double { return 1.0; };
  double boundary_functional =
      AvgValBoundary::compBoundaryFunctional(dofh, mu, w);

  std::cout << "H1s-norm: " << h1s_norm << "\n";
  std::cout << "F: " << boundary_functional << "\n";

#if SOLUTION
  auto results = AvgValBoundary::approxBoundaryFunctionalValues(7);
  double ground_truth = results[6].second;
  std::cout << "N_Dofs    Error\n";
  for (int i = 0; i < 6; i++) {
    std::cout << results[i].first << "  "
              << std::abs(results[i].second - ground_truth) << "\n";
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return 0;
}
/* SAM_LISTING_END_1 */
