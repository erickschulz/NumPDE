/**
 * @file regularizedneumannproblem_test.cc
 * @brief NPDE homework RegularizedNeumannProblem code
 * @author Christian Mitsch, Philippe Peter
 * @date March 2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <memory>
#include <utility>

// Lehrfem++ includes
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../getgalerkinlse.h"
#include "../regularizedneumannproblem.h"

namespace RegularizedNeumannProblem::test {

// Test for getGalerkinLSE method (not used in problem, but example code in PDF)
TEST(RegularizedNeumannProblem, getGalerkinLSE) {
  // We test on a triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto& dofh = fe_space_p->LocGlobMap();

  // Initialize constant mesh functions
  auto mf_f0 = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto mf_f2 = lf::mesh::utils::MeshFunctionConstant(2.0);
  auto mf_h0 = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto mf_h3 = lf::mesh::utils::MeshFunctionConstant(3.0);

  // Get boundary edges
  auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // Compute exact solution with f = 2.0, h=0.0
  Eigen::VectorXd f2_h0_solution(dofh.NumDofs());
  f2_h0_solution.setZero();

  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    auto geo_ptr = cell->Geometry();
    auto area = lf::geometry::Volume(*geo_ptr);
    auto glob_dof_indices = dofh.GlobalDofIndices(*cell);
    int local_dofs = dofh.NumLocalDofs(*cell);
    for (int i = 0; i < local_dofs; i++) {
      f2_h0_solution(glob_dof_indices[i]) += area / 3.0;
    }
  }

  f2_h0_solution *= 2.0;

  // Compute exact solution with f = 0.0, h=3.0
  Eigen::VectorXd f0_h3_solution(dofh.NumDofs());
  f0_h3_solution.setZero();

  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    if (bd_flags(*edge)) {
      auto geo_ptr = edge->Geometry();
      auto area = lf::geometry::Volume(*geo_ptr);
      auto glob_dof_indices = dofh.GlobalDofIndices(*edge);
      int local_dofs = dofh.NumLocalDofs(*edge);
      for (int i = 0; i < local_dofs; i++) {
        f0_h3_solution(glob_dof_indices[i]) += area / 2.0;
      }
    }
  }

  f0_h3_solution *= 3.0;

  // Get FEM Solution from getGalerkinLSE
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f2, mf_h0);
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair_2 =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f0, mf_h3);
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair_3 =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f2, mf_h3);

  // We only test rhs since A is simple laplacian
  auto approx_rhs_f2_h0 = sol_pair.second;
  auto approx_rhs_f0_h3 = sol_pair_2.second;
  auto approx_rhs_f2_h3 = sol_pair_3.second;

  EXPECT_NEAR((f2_h0_solution - approx_rhs_f2_h0).norm(), 0.0, 1.0e-8);
  EXPECT_NEAR((f0_h3_solution - approx_rhs_f0_h3).norm(), 0.0, 1.0e-8);

  // Test solution of f=2.0, h=3.0 -> Should be equal to adding up our two exact
  // solutions
  EXPECT_NEAR((f0_h3_solution + f2_h0_solution - approx_rhs_f2_h3).norm(), 0.0,
              1.0e-8);
}

// Test for sub-exercise c with functions h and f being constant
TEST(RegularizedNeumannProblem, solution_test_dropDof_const) {
  // read mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c =
      RegularizedNeumannProblem::getGalerkinLSE_dropDof(fe_space, f, h);

  // Now we compare the returned values to hard coded results
  const double eps = 1e-10;

  // Check Matrix
  Eigen::MatrixXd solution_mat_c(5, 5);
  solution_mat_c << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 1,
      -1, 0, -1, -1, -1, 4;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_c - result_c.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_c(5);
  solution_vec_c << 0, 1.1666666667, 1.1666666667, 1.1666666667, 0.33333333333;
  // Compare with expected results
  EXPECT_NEAR((solution_vec_c - result_c.second).norm(), 0.0, eps);
}

// Test for sub-exercise c with functions h and f not being constant
TEST(RegularizedNeumannProblem, solution_test_dropDof_gen) {
  // read  mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c =
      RegularizedNeumannProblem::getGalerkinLSE_dropDof(fe_space, f, h);

  // Now we compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_c(5, 5);
  solution_mat_c << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 1,
      -1, 0, -1, -1, -1, 4;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_c - result_c.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_c(5);
  solution_vec_c << 0, 1.1666666667, 1.91667, 1.1666666667, 0.33333333333;
  // Compare with expected results
  EXPECT_NEAR((solution_vec_c - result_c.second).norm(), 0.0, eps);
}

// Test for sub-exercise f with functions h and f being constant
TEST(RegularizedNeumannProblem, solution_test_augment_const) {
  // read mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_f =
      RegularizedNeumannProblem::getGalerkinLSE_augment(fe_space, f, h);

  // Compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_f(6, 6);
  solution_mat_f << 1, 0, 0, 0, -1, 0.1666667, 0, 1, 0, 0, -1, 0.1666667, 0, 0,
      1, 0, -1, 0.1666667, 0, 0, 0, 1, -1, 0.1666667, -1, -1, -1, -1, 4,
      0.3333333, 0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_f - result_f.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_f(6);
  solution_vec_f << 1.1666666667, 1.1666666667, 1.1666666667, 1.1666666667,
      0.33333333333, 0;
  EXPECT_NEAR((solution_vec_f - result_f.second).norm(), 0.0, eps);
}

// Test for sub-exercise f with functions h and f not being constant
TEST(RegularizedNeumannProblem, solution_test_augment_gen) {
  // read  mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_f =
      RegularizedNeumannProblem::getGalerkinLSE_augment(fe_space, f, h);

  // Compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_f(6, 6);
  solution_mat_f << 1, 0, 0, 0, -1, 0.1666667, 0, 1, 0, 0, -1, 0.1666667, 0, 0,
      1, 0, -1, 0.1666667, 0, 0, 0, 1, -1, 0.1666667, -1, -1, -1, -1, 4,
      0.3333333, 0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_f - result_f.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_f(6);
  solution_vec_f << 0.416667, 1.1666666667, 1.91667, 1.1666666667,
      0.33333333333, 0;
  EXPECT_NEAR((solution_vec_f - result_f.second).norm(), 0.0, eps);
}

}  // namespace RegularizedNeumannProblem::test
