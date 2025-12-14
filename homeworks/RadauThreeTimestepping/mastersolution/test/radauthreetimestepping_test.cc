/**
 * @file radauthreetimestepping_test.cc
 * @brief NPDE homework "RadauThreeTimestepping" code
 * @author Tobias Rohner, edited by Oliver Rietmann
 * @date 16.03.2020
 * @copyright Developed at ETH Zurich
 */

#include "../radauthreetimestepping.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace RadauThreeTimestepping::test {

TEST(RadauThreeTimestepping, TrapRuleLinFEElemVecProvider) {
  // Get some triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Define some easy functions to test the provider with
  auto f1 = [](const Eigen::Vector2d& x) { return 0.0; };
  auto f2 = [](const Eigen::Vector2d& x) { return 1.0; };
  auto f3 = [](const Eigen::Vector2d& x) { return x[0]; };
  // Check the element vector for each triangle
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f1p(f1);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f2p(f2);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f3p(f3);
  for (const auto tria : mesh_p->Entities(0)) {
    const auto geom = tria->Geometry();
    auto ev1 = f1p.Eval(*tria);
    auto ev2 = f2p.Eval(*tria);
    auto ev3 = f3p.Eval(*tria);
    ASSERT_TRUE(ev1.isApprox(Eigen::Vector3d::Zero()));
    ASSERT_TRUE(ev2.isApprox(
        Eigen::Vector3d::Constant(lf::geometry::Volume(*geom) / 3)));
    Eigen::Vector3d e3_correct =
        lf::geometry::Volume(*geom) / 3 * lf::geometry::Corners(*geom).row(0);
    ASSERT_TRUE(ev3.isApprox(e3_correct));
  }
}

TEST(RadauThreeTimestepping, rhsVectorheatSource) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto& dofh = fespace.LocGlobMap();
  // Assemble the vectors for t=0 and t=0.5
  const Eigen::VectorXd rhs0 =
      RadauThreeTimestepping::rhsVectorheatSource(dofh, 0.0);
  const Eigen::VectorXd rhs1 =
      RadauThreeTimestepping::rhsVectorheatSource(dofh, 0.5);
  // Get the DOFs on the boundary
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);

  // Create a functional for time t=0 and t=0.5
  auto f0 = [](const Eigen::Vector2d& x) {
    return ((x[0] - 0.5) * (x[0] - 0.5) + x[1] * x[1] < 0.25) ? 1.0 : 0.0;
  };
  auto f1 = [](const Eigen::Vector2d& x) {
    return (x[0] * x[0] + (x[1] - 0.5) * (x[1] - 0.5) < 0.25) ? 1.0 : 0.0;
  };

  // Assume TrapRuleLinFEElemVecProvider works correctly,
  // as it is tested above
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider provider0(f0);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider provider1(f1);
  Eigen::VectorXd rhs0_test = Eigen::VectorXd::Zero(dofh.NumDofs());
  Eigen::VectorXd rhs1_test = Eigen::VectorXd::Zero(dofh.NumDofs());
  lf::assemble::AssembleVectorLocally(0, dofh, provider0, rhs0_test);
  lf::assemble::AssembleVectorLocally(0, dofh, provider1, rhs1_test);

  double tol = 1e-10;

  // Make sure the TrapRuleLinFEElemVecProvider is indeed implemented already
  ASSERT_TRUE(rhs0_test.norm() > tol);
  ASSERT_TRUE(rhs1_test.norm() > tol);

  for (int i = 0; i < dofh.NumDofs(); ++i) {
    if (boundary(dofh.Entity(i))) {
      // Check whether boundary DOFs are set to zero
      ASSERT_DOUBLE_EQ(rhs0[i], 0)
          << "Have you forgotten to set the boundary values to zero?";
      ASSERT_DOUBLE_EQ(rhs1[i], 0)
          << "Have you forgotten to set the boundary values to zero?";
    } else {
      // Check whether the value coincides with the one
      // computed previously
      ASSERT_NEAR(rhs0[i], rhs0_test[i], tol);
      ASSERT_NEAR(rhs1[i], rhs1_test[i], tol);
    }
  }
}

TEST(RadauThreeTimestepping, solveHeatEvolution) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto& dofh = fespace.LocGlobMap();

  // Solve heat evolution with zero initial and boundary conditions
  double final_time = 1.0;
  unsigned int m = 50;

  Eigen::VectorXd sol =
      RadauThreeTimestepping::solveHeatEvolution(dofh, m, final_time);

  Eigen::VectorXd ref_sol(13);
  ref_sol << 0, 0, 0, 1.47965e-06, 1.46476e-06, 0, 0, 1.78839e-06, 1.36475e-06,
      0, 0, 0, 0;

  double tol = 1.e-4;

  ASSERT_NEAR((ref_sol - sol).norm(), 0.0, tol);
}

TEST(RadauThreeTimestepping, dropMatrixRowsColumns) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto& dofh = fespace.LocGlobMap();
  const lf::base::size_type N_dofs = dofh.NumDofs();

  // Obtain an array of boolean flags for the vertices of the mesh: 'true'
  // indicates that the vertex lies on the boundary.
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsColumns
  auto bdy_vertices_selector = [&bd_flags, &dofh](unsigned int idx) -> bool {
    return bd_flags(dofh.Entity(idx));
  };

  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix elLapMat_builder;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elLapMat_builder, A_COO);

  RadauThreeTimestepping::dropMatrixRowsColumns(bdy_vertices_selector, A_COO);

  Eigen::SparseMatrix<double> A_sps = A_COO.makeSparse();

  Eigen::MatrixXd A(A_sps);

  Eigen::MatrixXd A_ref = Eigen::MatrixXd::Zero(N_dofs, N_dofs);
  A_ref(0, 0) = 1;
  A_ref(1, 1) = 1;
  A_ref(2, 2) = 1;
  A_ref(3, 3) = 3.625;
  A_ref(3, 4) = -0.625;
  A_ref(3, 7) = -0.25;
  A_ref(3, 8) = -0.75;
  A_ref(4, 3) = -0.625;
  A_ref(4, 4) = 4.375;
  A_ref(4, 7) = -2;
  A_ref(5, 5) = 1;
  A_ref(6, 6) = 1;
  A_ref(7, 3) = -0.25;
  A_ref(7, 4) = -2;
  A_ref(7, 7) = 4.5;
  A_ref(7, 8) = -0.75;
  A_ref(8, 3) = -0.75;
  A_ref(8, 7) = -0.75;
  A_ref(8, 8) = 3.75;
  A_ref(9, 9) = 1;
  A_ref(10, 10) = 1;
  A_ref(11, 11) = 1;
  A_ref(12, 12) = 1;

  double tol = 1.e-4;

  ASSERT_NEAR((A - A_ref).norm(), 0.0, tol);
}

TEST(RadauThreeTimestepping, LinFEMassMatrixProvider) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto& dofh = fespace.LocGlobMap();

  // Compare the element matrices for every cell
  RadauThreeTimestepping::LinFEMassMatrixProvider provider;
  for (const auto cell : mesh_p->Entities(0)) {
    const auto geom = cell->Geometry();
    // Compute the correct element matrix
    Eigen::Matrix3d em_correct;
    em_correct.setConstant(lf::geometry::Volume(*geom) / 12);
    em_correct.diagonal() *= 2;
    // Get the element matrix from the element matrix provider
    Eigen::Matrix3d em_prov = provider.Eval(*cell);
    // Compare the matrices
    ASSERT_TRUE(em_correct.isApprox(em_prov));
  }
}

TEST(RadauThreeTimestepping, discreteEvolutionOperator) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto& dofh = fespace.LocGlobMap();

  // Compute A and M
  lf::uscalfe::LinearFELaplaceElementMatrix A_provider;
  RadauThreeTimestepping::LinFEMassMatrixProvider M_provider;
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::COOMatrix<double> M_COO(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A_provider, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_provider, M_COO);
  const auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);
  const auto selector = [&](unsigned idx) {
    return bd_flags(dofh.Entity(idx));
  };
  RadauThreeTimestepping::dropMatrixRowsColumns(selector, A_COO);
  RadauThreeTimestepping::dropMatrixRowsColumns(selector, M_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();

  const double dt = 0.0001;
  const Eigen::VectorXd phi =
      RadauThreeTimestepping::rhsVectorheatSource(dofh, dt);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(M + dt * A);
  RadauThreeTimestepping::Radau3MOLTimestepper timestepper(dofh);
  // Test for different mu
  for (unsigned i = 0; i < dofh.NumDofs(); ++i) {
    Eigen::VectorXd mu0(dofh.NumDofs());
    mu0.setZero();
    mu0[i] = 1;
    // Compute mu at the next timestep using discreteEvolutionOperator
    const Eigen::VectorXd mu_dEO =
        timestepper.discreteEvolutionOperator(0, dt, mu0);
    // Compute mu at the next timestep using implicit euler
    const Eigen::VectorXd mu_iE = solver.solve(M * mu0 + dt * phi);
    // Compare the two mu
    ASSERT_TRUE((mu_dEO - mu_iE).array().abs().maxCoeff() < 1e-4);
  }
}

}  // end namespace RadauThreeTimestepping::test
