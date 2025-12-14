/**
 * @file blendedparameterization_test.cc
 * @brief NPDE homework BlendedParameterization code
 * @author Bob Schreiner
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../blendedparameterization.h"

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace BlendedParameterization::test {
class BlendedParametrizationElementMatrixProvider {
 public:
  BlendedParametrizationElementMatrixProvider() = default;
  bool isActive(const lf::mesh::Entity&) { return true; }
  Eigen::MatrixXd Eval(const lf::mesh::Entity& cell) {
    LF_ASSERT_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                  "Only implemented for Triangles");
    auto geo_p = cell.Geometry();
    Eigen::MatrixXd coords = lf::geometry::Corners(*geo_p);
    const BlendedParameterization::coord_t a0 = coords.col(0);
    const BlendedParameterization::coord_t a1 = coords.col(1);
    const BlendedParameterization::coord_t a2 = coords.col(2);

    const BlendedParameterization::Segment gamma01(a0, a1);
    const BlendedParameterization::Segment gamma12(a1, a2);
    const BlendedParameterization::Segment gamma20(a2, a0);

    return BlendedParameterization::evalBlendLocMat(gamma01, gamma12, gamma20);
  }
};

TEST(BlendedParameterization, TestGalerkin) {
  // use test mesh (with only affine equivalent cells!) to set up fe space
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs(dofh.NumDofs());

  // Using the implemented function in blendedparametrization.cc to compute the
  // lhs Matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  BlendedParametrizationElementMatrixProvider elmat_builder;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // compute galerkin matrix using ReactionDiffusionElementMatrixProvider
  lf::assemble::COOMatrix<double> B(N_dofs, N_dofs);
  auto alpha = [](Eigen::Vector2d /*x*/) -> double { return 1.0; };
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
}

}  // namespace BlendedParameterization::test
