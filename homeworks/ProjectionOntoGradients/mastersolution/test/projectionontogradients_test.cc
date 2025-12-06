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
  // Building test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Divergence-free vector field
  const auto f = [](Eigen::Vector2d x) { return Eigen::Vector2d(-x(1), x(0)); };

  // DofHandler for the p.w. linear Lagrangian finite element
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1},
                                          {lf::base::RefEl::kSegment(), 0},
                                          {lf::base::RefEl::kTria(), 0}});
  // Compute solution
  const Eigen::VectorXd sol_vec =
      ProjectionOntoGradients::projectOntoGradients(dofh, f);

  // As stated in the exercise, we would expect the solution to be zero.
  // Hence we check every entry if its (numerically) zero.
  // The GoogleTest framework provides a function we may use:
  /*   EXPECT_NEAR(value 1, value 2, max. difference) */
  const double eps = 1e-15;
  for (std::size_t i = 0; i < sol_vec.size(); ++i) {
    EXPECT_NEAR(sol_vec[i], 0.0, eps);
    // Try testing for equality, you'll see it will fail miserably!
    /* EXPECT_EQ(sol_vec[i], 0.0); */
  }
}
/* SAM_LISTING_END_1 */

TEST(ProjectionOntoGradients, exact_sol_test) {
  // I. Construct the test mesh
  // mesh builder in a world of dimension 2
  lf::mesh::utils::TPTriagMeshBuilder my_builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));

  // define the test mesh
  my_builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(10)
      .setNumYCells(10);

  auto mesh_p = my_builder.Build();

  // II. Construct a linear finite element space on the test mesh
  lf::uscalfe::FeSpaceLagrangeO1<double> fe_space(mesh_p);

  // III. Define a function which computes the index of the triangle in which
  // the coorindates of a given point are
  const auto triangleIndex = [](Eigen::Vector2d x) {
    if (x(0) >= 0.0 && x(0) <= 0.5 && x(1) >= 0.0 && x(1) <= 0.5) {
      if (x(0) < x(1)) return 0;
      return 1;
    }
    if (x(0) >= 0.0 && x(0) <= 0.5 && x(1) >= 0.5 && x(1) <= 1.0) {
      if (x(0) + 0.5 < x(1)) return 2;
      return 3;
    }
    if (x(0) >= 0.5 && x(0) <= 1.0 && x(1) >= 0.0 && x(1) <= 0.5) {
      if (x(0) < 0.5 + x(1)) return 4;
      return 5;
    }
    if (x(0) >= 0.5 && x(0) <= 1.0 && x(1) >= 0.5 && x(1) <= 1.0) {
      if (x(0) < x(1)) return 6;
      return 7;
    }
    LF_ASSERT_MSG(false, "Coordinates outside of unit square");
    return -1;  // should never be reached, avoids warnings
  };

  // IV. Define the function f
  const auto f = [&triangleIndex](Eigen::Vector2d x) {
    int triang_idx = triangleIndex(x);
    // return the function value according to the definition in the
    // exercise
    switch (triang_idx) {
      case 0:
        return Eigen::Vector2d(2, 0);
      case 1:
        return Eigen::Vector2d(0, 2);
      case 3:
        return Eigen::Vector2d(2, -2);
      case 4:
        return Eigen::Vector2d(-2, 2);
      case 6:
        return Eigen::Vector2d(0, -2);
      case 7:
        return Eigen::Vector2d(-2, 0);
      default:
        return Eigen::Vector2d(0, 0);
    }
  };

  // V. Define the tent function associated with the central node
  const auto tentFunction = [&triangleIndex](Eigen::Vector2d x) {
    int triang_idx = triangleIndex(x);
    // Observe that the restriction of the tent function on
    // any adjecant triangle is by definition a linear function a*x + b*y + c.
    // The gradients of these linear functions are given in the exercise
    // so we can read of the values of a and b. We then compute the
    // value of c to ensure that the linear function evaluates to
    // one at the central node.
    switch (triang_idx) {
      case 0:
        return 2.0 * x(0);
      case 1:
        return 2.0 * x(1);
      case 3:
        return 2.0 * x(0) - 2.0 * x(1) + 1.0;
      case 4:
        return -2.0 * x(0) + 2.0 * x(1) + 1.0;
      case 6:
        return -2.0 * x(1) + 2.0;
      case 7:
        return -2.0 * x(0) + 2.0;
      default:
        return 0.0;
    }
  };

  // VI.Determine the coefficient vector of the tent function in
  // the FE space (perform a Nodal projection)
  // The Function lf::fe::NodalProjection requires a mesh function as its
  // second argument, so we first construct a meshFunction object which
  // describes the tent function.
  auto tentFunction_mf = lf::mesh::utils::MeshFunctionGlobal(tentFunction);
  auto ref_vec = lf::fe::NodalProjection<double>(fe_space, tentFunction_mf);

  // VII. Test your implementation
  const Eigen::VectorXd sol_vec =
      projectOntoGradients(fe_space.LocGlobMap(), f);
  EXPECT_NEAR((sol_vec - ref_vec).norm(), 0.0, 1e-12);
}

}  // namespace ProjectionOntoGradients::test
