/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author Erick Schulz, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include "../projectionontogradients.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <memory>

namespace ProjectionOntoGradients::test {

TEST(ProjectionOntoGradients, ElementMatrixProvider) {
  // Building triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // implemented element matrix provider
  ElementMatrixProvider my_elem_mat_provider{};

  // For comparison
  lf::uscalfe::LinearFELaplaceElementMatrix lfe_elem_mat_provider{};

  // loop over cells and compute element matrices
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    Eigen::Matrix3d my_mat{my_elem_mat_provider.Eval(*cell)};
    lf::uscalfe::LinearFELaplaceElementMatrix::ElemMat lfe_mat{
        lfe_elem_mat_provider.Eval(*cell)};

    // compare element matrices:
    EXPECT_NEAR((lfe_mat.block<3, 3>(0, 0) - my_mat).norm(), 0.0, 1E-3);
  }
}

TEST(ProjectionOntoGradients, GradProjRhsProvider_1) {
  // Building triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // initialize functions f and v for which the linear form is evaluated
  auto f = [](Eigen::Vector2d x) { return Eigen::Vector2d(1.0, 2.0); };
  auto v = [](Eigen::Vector2d x) { return 2 * x(0) + x(1); };
  auto v_mf = lf::mesh::utils::MeshFunctionGlobal(v);

  // set up finite elements
  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // construct element vector provider
  GradProjRhsProvider my_vec_provider(f);

  // assemble vector
  Eigen::VectorXd phi(fe_space->LocGlobMap().NumDofs());
  phi.setZero();
  lf::assemble::AssembleVectorLocally(0, fe_space->LocGlobMap(),
                                      my_vec_provider, phi);

  // project v onto the fe space
  auto v_vec = lf::fe::NodalProjection<double>(*fe_space, v_mf);

  // evaluate linear form on projected function:
  auto product = (v_vec.transpose() * phi).eval();
  EXPECT_NEAR(product(0, 0), 36, 1E-5);
}

TEST(ProjectionOntoGradients, GradProjRhsProvider_2) {
  // Building triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // initialize functions f and v for which the linear form is evaluated
  auto f = [](Eigen::Vector2d x) { return Eigen::Vector2d(x(1), x(0)); };
  auto v = [](Eigen::Vector2d x) { return 2 * x(0) + x(1); };
  auto v_mf = lf::mesh::utils::MeshFunctionGlobal(v);

  // set up finite elements
  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // construct element vector provider
  GradProjRhsProvider my_vec_provider(f);

  // assemble vector
  Eigen::VectorXd phi(fe_space->LocGlobMap().NumDofs());
  phi.setZero();
  lf::assemble::AssembleVectorLocally(0, fe_space->LocGlobMap(),
                                      my_vec_provider, phi);

  // project v onto the fe space
  auto v_vec = lf::fe::NodalProjection<double>(*fe_space, v_mf);

  // evaluate linear form on projected function:
  auto product = (v_vec.transpose() * phi).eval();
  EXPECT_NEAR(product(0, 0), 40.5, 1E-5);
}

/* SAM_LISTING_BEGIN_1 */
TEST(ProjectionOntoGradients, div_free_test) {
  //====================
  // Your code goes here
  //====================
}
/* SAM_LISTING_END_1 */

TEST(ProjectionOntoGradients, exact_sol_test) {
  // I. Construct the test mesh
  //====================
  // Your code goes here
  //====================

  // II. Construct a linear finite element space on the test mesh
  //====================
  // Your code goes here
  //====================

  // III. Define a function which computes the index of the triangle in which
  // the coorindates of a given point are
  //====================
  // Your code goes here
  //====================

  // IV. Define the function f
  //====================
  // Your code goes here
  //====================

  // V. Define the tent function associated with the central node
  //====================
  // Your code goes here
  //====================

  // VI.Determine the coefficient vector of the tent function in
  // the FE space (perform a Nodal projection)
  //====================
  // Your code goes here
  //====================

  // VII. Test your implementation
  //====================
  // Your code goes here
  //====================
}

}  // namespace ProjectionOntoGradients::test
