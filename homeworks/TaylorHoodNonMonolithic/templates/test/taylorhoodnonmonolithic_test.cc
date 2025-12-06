/**
 * @file TaylorHoodNonMonolithic_test.cc
 * @brief NPDE homework TaylorHoodNonMonolithic code
 * @author Ralf Hiptmair
 * @date July 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../taylorhoodnonmonolithic.h"

#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <gtest/gtest.h>
#include <lf/base/ref_el.h>
#include <lf/fe/fe_tools.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>

#include <Eigen/Core>
#include <cstddef>

/* Test in the google testing framework

  The following assertions are available, syntax
  EXPECT_XX( ....) << [anything that can be givne to std::cerr]

  EXPECT_EQ(val1, val2)
  EXPECT_NEAR(val1, val2, abs_error) -> should be used for numerical results!
  EXPECT_NE(val1, val2)
  EXPECT_TRUE(condition)
  EXPECT_FALSE(condition)
  EXPECT_GE(val1, val2)
  EXPECT_LE(val1, val2)
  EXPECT_GT(val1, val2)
  EXPECT_LT(val1, val2)
  EXPECT_STREQ(str1,str2)
  EXPECT_STRNE(str1,str2)
  EXPECT_STRCASEEQ(str1,str2)
  EXPECT_STRCASENE(str1,str2)

  "EXPECT" can be replaced with "ASSERT" when you want to program to terminate,
 if the assertion is violated.
 */

namespace TaylorHoodNonMonolithic::test {

TEST(TaylorHoodNonMonolithic, THBElementMatrixProvider1) {
  // Define a pre-computed element matrix
  Eigen::MatrixXd exact_elem(3, 6);
  exact_elem << -0.0555556, 0, 0, 0.0555556, 0.0555556, -0.0555556, 0,
      0.0555556, 0, -0.0555556, 0.0555556, -0.0555556, 0, 0, 0, 0, 0.111111,
      -0.111111;

  // Generate the mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto cells = mesh_ptr->Entities(0);

  // Compute the element matrix for the first cell
  TaylorHoodNonMonolithic::THBElementMatrixProvider THBemp(X_Dir);
  auto elem = THBemp.Eval(*cells[0]);

  // Check if the solution is correct
  ASSERT_EQ(elem.cols(), 6);
  ASSERT_EQ(elem.rows(), 3);
  ASSERT_NEAR((elem - exact_elem).norm(), 0., 1e-4);
}

TEST(TaylorHoodNonMonolithic, THBElementMatrixProvider2) {
  // Define a pre-computed element matrix
  Eigen::MatrixXd exact_elem(3, 6);
  exact_elem << -0.0277778, 0, 0, -0.138889, 0.0277778, 0.138889, 0, -0.0555556,
      0, -0.111111, 0.111111, 0.0555556, 0, 0, 0.0833333, -0.0833333,
      -0.0277778, 0.0277778;

  // Generate the mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto cells = mesh_ptr->Entities(0);

  // Compute the element matrix for the first cell
  TaylorHoodNonMonolithic::THBElementMatrixProvider THBemp(Y_Dir);
  auto elem = THBemp.Eval(*cells[0]);

  // Check if the solution is correct
  ASSERT_EQ(elem.cols(), 6);
  ASSERT_EQ(elem.rows(), 3);
  ASSERT_NEAR((elem - exact_elem).norm(), 0., 1e-4);
}

TEST(TaylorHoodNonMonolithic, Uzawa) {
  // Obtain mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Generate Galerkin matrices
  const size_t N = 10;
  const size_t M = 5;
  Eigen::SparseMatrix<double> A(N, N);
  for (int k = 0; k < N; ++k) {
    A.insert(k, k) = 2.0;
  }
  for (int k = 1; k < N; ++k) {
    A.insert(k - 1, k) = -1.0;
    A.insert(k, k - 1) = -1.0;
  }
  Eigen::SparseMatrix<double> B_x(M, N);
  Eigen::SparseMatrix<double> B_y(M, N);
  for (int k = 0; k < M; ++k) {
    B_x.insert(k, k) = 1.0;
    B_y.insert(k, N - 1 - k) = 1.0;
  }
  std::cout << "A = \n"
            << Eigen::MatrixXd(A) << "\n B_x = \n"
            << Eigen::MatrixXd(B_x) << "\n B_y = \n"
            << Eigen::MatrixXd(B_y) << std::endl;

  Eigen::VectorXd phi_x = Eigen::VectorXd::Constant(N, 1.0);
  Eigen::VectorXd phi_y = Eigen::VectorXd::LinSpaced(N, 0.0, 1.0);
  // For recording progress of the iteration
  std::vector<std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>>
      rec_data;
  auto rec = [&rec_data](const Eigen::VectorXd& mu_x,
                         const Eigen::VectorXd& mu_y,
                         const Eigen::VectorXd& pi) -> void {
    rec_data.emplace_back(mu_x, mu_y, pi);
  };
  [[maybe_unused]] auto [res_mu_x, res_mu_y, res_pi] =
      TaylorHoodNonMonolithic::CGUzawa(A, B_x, B_y, phi_x, phi_y, 1E-6, 1E-8,
                                       10, rec);
  std::cout << "CG Uzawa took " << rec_data.size() << "steps\n";
  int step = 1;
  for (auto vecs : rec_data) {
    auto& [mu_x, mu_y, pi] = vecs;
    const Eigen::VectorXd res_x = phi_x - A * mu_x - B_x.transpose() * pi;
    const Eigen::VectorXd res_y = phi_y - A * mu_y - B_y.transpose() * pi;
    const Eigen::VectorXd res_pi = B_x * mu_x + B_y * mu_y;
    EXPECT_NEAR(res_x.norm(), 0.0, 1.0E-10);
    EXPECT_NEAR(res_y.norm(), 0.0, 1.0E-10);
    std::cout << "step " << step << ": |res_x| = " << res_x.norm()
              << ", |res_y| = " << res_y.norm()
              << ", |res_pi| = " << res_pi.norm() << std::endl;
    if (step == 10) {
      Eigen::VectorXd exact_final_mu_x(10), exact_final_mu_y(10);
      exact_final_mu_x << 0.833333, 1.66667, 2.44444, 3.11111, 3.61111, 5.50926,
          6.40741, 6.30556, 5.2037, 3.10185;
      exact_final_mu_y << -0.231481, -0.462963, -0.805556, -1.37037, -2.26852,
          -3.61111, -3.11111, -2.44444, -1.66667, -0.833333;
      EXPECT_NEAR((mu_x - exact_final_mu_x).norm(), 0., 1e-6);
      EXPECT_NEAR((mu_y - exact_final_mu_y).norm(), 0., 1e-6);
    }
    step++;
  }
}

}  // namespace TaylorHoodNonMonolithic::test
