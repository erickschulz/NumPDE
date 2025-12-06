/**
 * @file zienkiewiczzhuestimator_test.cc
 * @brief NPDE homework "ZienkiewiczZhuEstimator" code
 * @author Philipp Lindenberger
 * @date 25.03.2020
 * @copyright Developed at ETH Zurich
 */

#include <memory>
// Eigen includes
#include <gtest/gtest.h>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <unsupported/Eigen/KroneckerProduct>

#include "../zienkiewiczzhuestimator.h"

namespace ZienkiewiczZhuEstimator::test {

/**
 * @brief test VectorProjectionMatrixProvider implementation
 */
TEST(ZienkiewiczZhuEstimator, VectorProjectionMatrixProvider) {
  // Load triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Pointer to scalar FE-space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::Matrix2d identity = Eigen::MatrixXd::Identity(2, 2);
  auto mf_one = lf::mesh::utils::MeshFunctionConstant(1.0);
  auto mf_zero = lf::mesh::utils::MeshFunctionConstant(0.0);

  // Student solution
  ZienkiewiczZhuEstimator::VectorProjectionMatrixProvider student_elmat_builder;
  // Reference solution
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_zero),
                                                      decltype(mf_one)>
      elmat_builder_exact(fe_space_p, mf_zero, mf_one);

  // Loop over all dim-2 entities, compute error
  for (const auto tria : mesh_p->Entities(0)) {
    auto student_elem_mat = student_elmat_builder.Eval(*tria);
    auto elem_mat_scalar_exact = elmat_builder_exact.Eval(*tria);
    auto elem_mat_exact =
        Eigen::KroneckerProduct(elem_mat_scalar_exact, identity);
    double error = (student_elem_mat - elem_mat_exact).norm();
    EXPECT_LT(error, 1.0e-12);
  }
}

/**
 * @brief test GradientProjectionVectorProvider implementation
 */
TEST(ZienkiewiczZhuEstimator, GradientProjectionVectorProvider) {
  // Load triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Pointer to scalar FE-space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Scalar DOF-Handler
  auto& dofh = fe_space_p->LocGlobMap();

  // Solution vector 1:
  auto mu_x = [](Eigen::Vector2d x) -> double { return x[0]; };
  auto mf_mu_x = lf::mesh::utils::MeshFunctionGlobal(
      [&mu_x](Eigen::Vector2d x) -> double { return mu_x(x); });
  auto mu_x_vec = lf::fe::NodalProjection(*fe_space_p, mf_mu_x);

  // Solution vector 2:
  auto mu_y = [](Eigen::Vector2d x) -> double { return x[1]; };
  auto mf_mu_y = lf::mesh::utils::MeshFunctionGlobal(
      [&mu_y](Eigen::Vector2d x) -> double { return mu_y(x); });
  auto mu_y_vec = lf::fe::NodalProjection(*fe_space_p, mf_mu_y);

  // Retrieve unit triangle from mesh (see documentation of
  // GenerateHybrid2DTestMesh(3) on LF++)
  auto unit_tria = mesh_p->EntityByIndex(0, 1);

  // Student solutions
  ZienkiewiczZhuEstimator::GradientProjectionVectorProvider grad_vec_builder_x(
      fe_space_p, mu_x_vec);
  ZienkiewiczZhuEstimator::GradientProjectionVectorProvider grad_vec_builder_y(
      fe_space_p, mu_y_vec);

  auto unit_tria_grad_x = grad_vec_builder_x.Eval(*unit_tria);
  auto unit_tria_grad_y = grad_vec_builder_y.Eval(*unit_tria);

  for (int i = 0; i < 6; i += 2) {
    ASSERT_NEAR(unit_tria_grad_x(i), 0.25, 1.0e-8);
    ASSERT_NEAR(unit_tria_grad_x(i + 1), 0.0, 1.0e-8);
    ASSERT_NEAR(unit_tria_grad_y(i), 0.0, 1.0e-8);
    ASSERT_NEAR(unit_tria_grad_y(i + 1), 0.25, 1.0e-8);
  }
}

/**
 * @brief test computeLumpedProjection implementation
 */
TEST(ZienkiewiczZhuEstimator, computeLumpedProjection) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain reference to scalar dofh
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Produce a dof handler for the vector-valued finite element space
  lf::assemble::UniformFEDofHandler vec_dofh(mesh_p,
                                             {{lf::base::RefEl::kPoint(), 2},
                                              {lf::base::RefEl::kSegment(), 0},
                                              {lf::base::RefEl::kTria(), 0},
                                              {lf::base::RefEl::kQuad(), 0}});

  // Solution vector 1:
  auto mu_x = [](Eigen::Vector2d x) -> double { return x[0]; };
  auto mf_mu_x = lf::mesh::utils::MeshFunctionGlobal(
      [&mu_x](Eigen::Vector2d x) -> double { return mu_x(x); });
  auto mu_x_vec = lf::fe::NodalProjection(*fe_space_p, mf_mu_x);

  // Solution vector 2:
  auto mu_y = [](Eigen::Vector2d x) -> double { return x[1]; };
  auto mf_mu_y = lf::mesh::utils::MeshFunctionGlobal(
      [&mu_y](Eigen::Vector2d x) -> double { return mu_y(x); });
  auto mu_y_vec = lf::fe::NodalProjection(*fe_space_p, mf_mu_y);

  // Student solutions
  auto sol_x = ZienkiewiczZhuEstimator::computeLumpedProjection(dofh, mu_x_vec,
                                                                vec_dofh);
  auto sol_y = ZienkiewiczZhuEstimator::computeLumpedProjection(dofh, mu_y_vec,
                                                                vec_dofh);

  for (int i = 0; i < sol_x.size(); i += 2) {
    ASSERT_NEAR(sol_x(i), 1.0, 1.0e-8);
    ASSERT_NEAR(sol_x(i + 1), 0.0, 1.0e-8);
    ASSERT_NEAR(sol_y(i), 0.0, 1.0e-8);
    ASSERT_NEAR(sol_y(i + 1), 1.0, 1.0e-8);
  }
}

/**
 * @brief test computeL2Deviation implementation
 */
TEST(ZienkiewiczZhuEstimator, computeL2Deviation) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain reference to scalar dofh
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Produce a dof handler for the vector-valued finite element space
  lf::assemble::UniformFEDofHandler vec_dofh(mesh_p,
                                             {{lf::base::RefEl::kPoint(), 2},
                                              {lf::base::RefEl::kSegment(), 0},
                                              {lf::base::RefEl::kTria(), 0},
                                              {lf::base::RefEl::kQuad(), 0}});

  auto mu = ZienkiewiczZhuEstimator::solveBVP(fe_space_p);

  // Assumes computeLumpedProjection works as expected (tested above)
  auto mu_grad =
      ZienkiewiczZhuEstimator::computeLumpedProjection(dofh, mu, vec_dofh);
  // Student solution
  double error_grad = computeL2Deviation(dofh, mu, vec_dofh, mu_grad);

  ASSERT_NEAR(error_grad, 0.0419009, 1.0e-6);
}

}  // namespace ZienkiewiczZhuEstimator::test
