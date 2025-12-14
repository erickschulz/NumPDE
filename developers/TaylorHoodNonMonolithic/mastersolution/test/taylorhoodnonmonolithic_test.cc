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
#if SOLUTION
TEST(TaylorHoodNonMonolithic, MatrixTests) {
  // Obtain a triangular mesh of the unit square
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(10)
      .setNumYCells(10);
  auto mesh_p = builder.Build();
  // Force functor (rotational force)
  auto force = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return {-x[1], x[0]};
  };
  // Generate Galerkin matrices
  auto [A, B_x, B_y, phi_x, phi_y] =
      TaylorHoodNonMonolithic::buildStokesLSE(mesh_p, force);

  // Quadratic Lagrange FEM:
  // This test takes for granted that this object uses the same loca -> global
  // map the object of the same type used inside buildStokesLSE()
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh_LO2{fe_space->LocGlobMap()};
  const lf::assemble::size_type n_LO2 = dofh_LO2.NumDofs();
  auto fe_LO1 =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh_LO1{fe_LO1->LocGlobMap()};
  const lf::assemble::size_type n_LO1 = dofh_LO1.NumDofs();
  {
    std::cout << "Testing matrix A\n";
    // A harmonic quadratic function
    auto harmquad = [](Eigen::Vector2d x) -> double {
      return (x[0] * x[0] - x[1] * x[1]);
    };
    // Create FE LO2 representation
    lf::mesh::utils::MeshFunctionGlobal hq_mf(harmquad);
    const Eigen::VectorXd hq_vec{lf::fe::NodalProjection(*fe_space, hq_mf)};
    // Multiply with matrix A
    const Eigen::VectorXd res_vec = A * hq_vec;
    // The entries of the result vector correspomnding to d.o.f.s in the
    // interior should vanish
    for (lf::assemble::glb_idx_t dof_idx = 0; dof_idx < n_LO2; ++dof_idx) {
      const lf::mesh::Entity& ent = dofh_LO2.Entity(dof_idx);
      const lf::geometry::Geometry& geo = *ent.Geometry();
      const Eigen::MatrixXd corners{lf::geometry::Corners(geo)};
      const Eigen::VectorXd mp{corners.rowwise().sum() / corners.cols()};
      if ((mp[0] > 1.0 / 3.0) and (mp[0] < 2.0 / 3.0) and
          (mp[1] > 1.0 / 3.0) and (mp[1] < 2.0 / 3.0)) {
        EXPECT_NEAR(res_vec[dof_idx], 0.0, 1.0E-6)
            << ", LO2 dof_idx = " << dof_idx << "\n";
      }
    }
  }
  {
    std::cout << "Testing B matrices\n";
    // Components of a divergence-free quadratic vector field
    lf::mesh::utils::MeshFunctionGlobal solen_x{
        [](Eigen::Vector2d x) -> double {
          return (-x[0] * x[0] + 3 * x[1] * x[1]);
        }};

    lf::mesh::utils::MeshFunctionGlobal solen_y{
        [](Eigen::Vector2d x) -> double {
          return (3 * x[0] * x[0] + 2 * x[0] * x[1]);
        }};
    const Eigen::VectorXd solen_x_vec{
        lf::fe::NodalProjection(*fe_space, solen_x)};
    const Eigen::VectorXd solen_y_vec{
        lf::fe::NodalProjection(*fe_space, solen_y)};
    const Eigen::VectorXd res_vec = B_x * solen_x_vec + B_y * solen_y_vec;
    // "Interior" components of result vector should vanish
    for (lf::assemble::glb_idx_t dof_idx = 0; dof_idx < n_LO1; ++dof_idx) {
      const lf::mesh::Entity& ent = dofh_LO1.Entity(dof_idx);
      EXPECT_TRUE(ent.RefEl() == lf::base::RefEl::kPoint());
      const lf::geometry::Geometry& geo = *ent.Geometry();
      const Eigen::MatrixXd corners{lf::geometry::Corners(geo)};
      const Eigen::VectorXd mp{corners.rowwise().sum() / corners.cols()};
      if ((mp[0] > 1.0 / 3.0) and (mp[0] < 2.0 / 3.0) and
          (mp[1] > 1.0 / 3.0) and (mp[1] < 2.0 / 3.0)) {
        EXPECT_NEAR(res_vec[dof_idx], 0.0, 1.0E-6)
            << ", LO1 dof_idx = " << dof_idx << "\n";
      }
    }
  }
  {
    std::cout << "Testing kernel of transposed B matrices\n";
    EXPECT_NEAR(
        (B_x.transpose() * Eigen::VectorXd::Constant(n_LO1, 1.0)).norm(), 0.0,
        1.0E-6);
    EXPECT_NEAR(
        (B_y.transpose() * Eigen::VectorXd::Constant(n_LO1, 1.0)).norm(), 0.0,
        1.0E-6);
  }
}
#endif

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
