/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author Erick Schulz, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

#include "projectionontogradients.h"

int main() {
  // for this exercise a main file is not required
  // but feel free to use it to call some of your functions for debugging
  // purposes
  const auto f = [](Eigen::Vector2d x) { return Eigen::Vector2d(-x(1), x(0)); };
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // Compute solution
  const Eigen::VectorXd sol_vec =
      ProjectionOntoGradients::projectOntoGradients(dofh, f);
  std::cout << sol_vec << std::endl;

  std::cout << "You may use this main file to call your function "
               "ProjectionOntoGradients::projectOntoGradients"
            << std::endl;
}
