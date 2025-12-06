/**
 * @ file PointEvaluationRhs_test.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author ?, Liaowang Huang (refactoring)
 * @ date ?, 06/01/2020 (refactoring)
 * @ copyright Developed at ETH Zurich
 */

#include "../pointevaluationrhs.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>
#include <iostream>

#include "../pointevaluationrhs_norms.h"

/* SAM_LISTING_BEGIN_1 */
void testGlobalInverseQuad(const lf::mesh::Entity& quad, Eigen::Vector2d xh) {
  LF_ASSERT_MSG(quad.RefEl() == lf::base::RefEl::kQuad(),
                "Cell must be a quadrilateral");
  // get the coordinates of the corners of this cell
  const lf::geometry::Geometry* geo_ptr = quad.Geometry();
  auto vertices = lf::geometry::Corners(*geo_ptr);
  // Image of point in unit square under parametric mapping
  Eigen::Vector2d x = geo_ptr->Global(xh);
  Eigen::Vector2d xh_comp = PointEvaluationRhs::GlobalInverseQuad(vertices, x);
  EXPECT_NEAR((xh - xh_comp).norm(), 0.0, 1.0E-8)
      << "quadl " << quad << ": Mismatch xh = " << xh << ", x = " << x
      << ", xh_comp = " << xh_comp << std::endl;
}
/* SAM_LISTING_END_1 */

// This test checks whether the implemented inverse mappings
// of the transformation from the reference element work correctly
// for the cells of a general hybrid mesh.
TEST(PoinEvaluationRhs, mapping_test) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Point in reference element for testing
  Eigen::Vector2d xh(0.27, 0.41);

  for (auto cell : mesh_p->Entities(0)) {
    // Get shape of cell
    const lf::geometry::Geometry* geo_ptr = cell->Geometry();
    // Get cordinates of vertices
    auto vertices = lf::geometry::Corners(*geo_ptr);
    // Global coordinates of testing point
    Eigen::Vector2d x{geo_ptr->Global(xh)};
    // Reference coordinates
    Eigen::Vector2d xh_comp;
    // Query type of cell
    const lf::base::RefEl ref_el = cell->RefEl();
    // Depending on the type of cell compute the pre-image of
    // the point x
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        xh_comp = PointEvaluationRhs::GlobalInverseTria(vertices, x);
        break;
      }
      case lf::base::RefEl::kQuad(): {
        xh_comp = PointEvaluationRhs::GlobalInverseQuad(vertices, x);
        break;
      }
      default: {
        LF_ASSERT_MSG(false, "Not implemented for " << ref_el);
        break;
      }
    }  // end switch

    ASSERT_NEAR((xh - xh_comp).norm(), 0.0, 1.0E-8);
  }  // end loop over cells
}

TEST(PoinEvaluationRhs, solution_test) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

  Eigen::VectorXd sol_vec;
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1}});
  auto result = PointEvaluationRhs::normsSolutionPointLoadDirichletBVP(
      dofh, Eigen::Vector2d(1.3, 1.7), sol_vec);

  double eps = 1e-6;

  // Test norms
  EXPECT_NEAR(result.first, 0.290660239, eps);
  EXPECT_NEAR(result.second, 0.39855911, eps);

  // Test solution
  Eigen::VectorXd correct_sol(15);

  correct_sol << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.166856, 0.0304501, 0.0521068,
      0.10716, 0.172936, 0.152714;

  for (int i = 0; i < sol_vec.size(); ++i) {
    EXPECT_NEAR(sol_vec(i), correct_sol(i), eps);
  }
}
