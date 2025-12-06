/**
 * @file
 * @brief UNITESTS for NPDE homework OutputImpedanceBVP
 * @author Erick Schulz
 * @date 28/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "../outputimpedancebvp.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <string>

#include "../evalclass.h"

namespace OutputImpedanceBVP::test {

TEST(OutputImpedanceBVP, computeApproxSolDirichlet) {
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/unitsquare.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};

  // Exact solution and Dirichlet boundary conditions
  Eigen::Vector2d g;
  g << 1.0, 3.0;
  auto uExact = [&g](Eigen::Vector2d x) -> double { return g.dot(x); };
  auto uExact_vec = interpolateData<std::function<double(Eigen::Vector2d)>>(
      fe_space_p, std::move(uExact));

  // Solve BVP
  Eigen::VectorXd uApprox_vec = solveImpedanceBVP(fe_space_p, g);

  ASSERT_TRUE(uApprox_vec.isApprox(uExact_vec));
}

}  // namespace OutputImpedanceBVP::test
