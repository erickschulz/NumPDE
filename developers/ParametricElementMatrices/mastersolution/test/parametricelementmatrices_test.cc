#include <gtest/gtest.h>
#include <lf/fe/fe.h>

#include "../anisotropicdiffusionelementmatrixprovider.h"
#include "../fesourceelemvecprovider.h"

namespace ParametricElementMatrices::test {

/* SAM_LISTING_BEGIN_1 */
TEST(ParametricElementMatrices, TestGalerkin) {
#if SOLUTION
  // use test mesh (with only affine equivalent cells!) to set up fe space
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(5, 1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs(dofh.NumDofs());

  // compute galerkin matrix for d(x) = sin(|x|)x using the implemented class
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  auto d = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (x * std::sin(x.norm()));
  };
  auto elmat_builder =
      ParametricElementMatrices::AnisotropicDiffusionElementMatrixProvider(d);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  // compute galerkin matrix using ReactionDiffusionElementMatrixProvider
  lf::assemble::COOMatrix<double> B(N_dofs, N_dofs);
  auto alpha = [&d](Eigen::Vector2d x) -> Eigen::Matrix2d {
    const Eigen::Vector2d d_x{d(x)};
    return Eigen::Matrix2d::Identity(2, 2) + d_x * d_x.transpose();
  };
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
  auto zero = [](Eigen::Vector2d x) -> double { return 0.; };
  lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};
  // set up quadrature rule to be able to compare
  std::map<lf::base::RefEl, lf::quad::QuadRule> quad_rules{
      {lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()},
      {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_EdgeMidpointRule()}};

  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_alpha), decltype(mf_zero)>
      elmat_builder_org(fe_space, mf_alpha, mf_zero, quad_rules);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder_org, B);

  auto A_crs = A.makeSparse();
  auto B_crs = B.makeSparse();
  // compare results (floating point comparison!)
  for (int i = 0; i < N_dofs; i++) {
    for (int j = 0; j < N_dofs; j++) {
      ASSERT_NEAR(A_crs.coeff(i, j), B_crs.coeff(i, j), 1E-9);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
TEST(ParametricElementMatrices, TestLoad) {
#if SOLUTION
  // use test mesh (with only affine equivalent cells!) to set up fe space
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(5, 1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs(dofh.NumDofs());
  // An affine linear function that can be represented exactly
  // in the space of p.w. linear Lagrangian finite element functions
  auto w_func = [](Eigen::Vector2d x) -> double { return (3.0 * x[0] - x[1]); };
  lf::mesh::utils::MeshFunctionGlobal mf_w_func{w_func};
  // function f(x) = 1/(1 + w(x)^2)
  auto f = [&w_func](Eigen::Vector2d x) -> double {
    return 1. / (1. + std::pow(w_func(x), 2));
  };
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // interpolation of w on fe space: reproduces function
  auto w = lf::fe::NodalProjection<double>(*fe_space, mf_w_func);

  // use own method to assemble load vector
  auto elvec_builder =
      ParametricElementMatrices::FESourceElemVecProvider(fe_space, w);
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  AssembleVectorLocally(0, dofh, elvec_builder, phi);
  // set up midpoint quadrature rules for all cell types
  std::map<lf::base::RefEl, lf::quad::QuadRule> quad_rules{
      {lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()},
      {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_EdgeMidpointRule()}};

  // compute library solution
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_lib(N_dofs);
  phi_lib.setZero();
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder_(fe_space, mf_f, quad_rules);
  AssembleVectorLocally(0, dofh, elvec_builder_, phi_lib);

  // compare the results
  for (int i = 0; i < N_dofs; i++) {
    ASSERT_NEAR(phi(i), phi_lib(i), 1E-9);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
}
/* SAM_LISTING_END_2 */

}  // namespace ParametricElementMatrices::test
