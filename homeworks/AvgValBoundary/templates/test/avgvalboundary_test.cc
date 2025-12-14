/**
 * @ file avgvalboundary_test.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans, edited by Oliver Rietmann
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "../avgvalboundary.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <memory>

namespace AvgValBoundary::test {

constexpr char mesh_file[] = "meshes/square.msh";

constexpr auto const_one = [](Eigen::Vector2d x) -> double { return 1.0; };

TEST(AvgValBoundary, TestH1SemiNorm) {
  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = solveTestProblem(dofh);
  // compute H1 seminorm
  double h1s_norm = compH1seminorm(dofh, mu);

  ASSERT_NEAR(h1s_norm, 0.151178, 0.00001);
}

TEST(AvgValBoundary, TestBoundaryFunctional) {
  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = solveTestProblem(dofh);
  // compute boundary functional
  double boundary_functional = compBoundaryFunctional(dofh, mu, const_one);
  ASSERT_NEAR(boundary_functional, 0.77546, 0.00001);
}
}  // namespace AvgValBoundary::test
