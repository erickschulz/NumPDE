/**
 * @file
 * @brief NPDE homework TEMPLATE MAIN FILE
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../hierarchicalerrorestimator.h"

#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <iostream>

namespace HEST::test {

TEST(HEST, trfLinToQuad) {
  /* Macros available in the Google test framework:
     EXPECT_EQ(x,y), EXPECT_NE(x,y), EXPECT_LT(x,y), EXPECT_LE(x,y),
     EXPECT_GT(x,y), EXPECT_GE(x,y) EXPECT_STREQ(x,y), EXPECT_STRNE(x,y) -> for
     C-strings only ! EXPECT_NEAR(x,y,abs_tol) All testing macros can output a
     message by a trailing << ....
   */
  // Obtain test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO2<double>> quad_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_ptr);
  const lf::uscalfe::FeSpaceLagrangeO2<double>& quad_space{*quad_space_p};
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> lfe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_ptr);
  const lf::uscalfe::FeSpaceLagrangeO1<double>& lfe_space{*lfe_space_p};
  // Get references to DofHandlers
  const lf::assemble::DofHandler& dh_quad{quad_space.LocGlobMap()};
  const lf::assemble::DofHandler& dh_lfe{lfe_space.LocGlobMap()};
  // Set up coefficient vector
  const lf::base::size_type N_dofs(dh_lfe.NumDofs());
  Eigen::VectorXd mu{Eigen::VectorXd::LinSpaced(N_dofs, 0.0, 1.0)};
  // Call conversion function
  Eigen::VectorXd nu{trfLinToQuad(lfe_space_p, quad_space_p, mu)};
  // Check agreement of L2 norm
  const lf::fe::MeshFunctionFE mf_quad(quad_space_p, nu);
  const lf::fe::MeshFunctionFE mf_lfe(lfe_space_p, mu);
  double L2diff =  // NOLINT
      std::sqrt(lf::fe::IntegrateMeshFunction(
          *mesh_ptr, lf::mesh::utils::squaredNorm(mf_quad - mf_lfe), 4));
  std::cout << "TEST trfLinToQuad: L2-norm of difference = " << L2diff
            << std::endl;
  EXPECT_NEAR(L2diff, 0.0, 1E-6);
}

TEST(HEST, compHierSurplusSolution) {
  // Obtain test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  const lf::mesh::Mesh& mesh{*mesh_ptr};
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO2<double>> quad_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_ptr);
  const lf::uscalfe::FeSpaceLagrangeO2<double>& quad_space{*quad_space_p};
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> lfe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_ptr);
  const lf::uscalfe::FeSpaceLagrangeO1<double>& lfe_space{*lfe_space_p};
  // Get references to DofHandlers
  const lf::assemble::DofHandler& dh_quad{quad_space.LocGlobMap()};
  const lf::assemble::DofHandler& dh_lfe{lfe_space.LocGlobMap()};
  // Set up coefficient vector
  const lf::base::size_type N_dofs(dh_lfe.NumDofs());
  Eigen::VectorXd mu{Eigen::VectorXd::LinSpaced(N_dofs, 0.0, 1.0)};
  // Diffusion coefficient and right-hand side as MeshFunctions
  auto alpha = [](Eigen::Vector2d /*x*/) -> Eigen::Matrix<double, 2, 2> {
    return Eigen::Matrix<double, 2, 2>::Identity();
  };
  auto f = [](Eigen::Vector2d /*x*/) -> double { return -4.0; };
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};
  // Do the real thing
  Eigen::VectorXd surplus =
      compHierSurplusSolution(mf_alpha, mf_f, lfe_space_p, quad_space_p, mu);
  // The first #nodes components of the result vector should vanish
  for (int k = 0; k < dh_lfe.NumDofs(); ++k) {
    EXPECT_NEAR(surplus[k], 0.0, 1E-10);
  }
  // Still no idea how to test the other components
}

}  // namespace HEST::test
