/**
 * @file MagStat2D_test.cc
 * @brief NPDE homework MagStat2D code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../magstat2d.h"

#include <gtest/gtest.h>
#include <lf/fe/fe_tools.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/mesh_function_unary.h>
#include <lf/uscalfe/fe_space_lagrange_o2.h>

#include <Eigen/Core>

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

namespace MagStat2D::test {

TEST(MagStat2D, PleaseNameTest) {
  int refsteps = 4;
  using namespace std::numbers;
  // Constant coefficient mu = 1
  lf::mesh::utils::MeshFunctionConstant<double> mf_one(1.0);
  // ********** Part I: Manufactured solution **********
  // Domain: unit square
  // A divergence-free vector field with zero tangential components
  // on the boundary of unit sqaure
  auto u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-std::cos(pi * x[0]) * std::sin(pi * x[1]),
                           std::sin(pi * x[0]) * std::cos(pi * x[1]));
  };
  // The solution component p is zero in this case
  auto p = [](Eigen::Vector2d /*x*/) -> double { return 0.0; };

  // Right-hand side source current; u is an eigenfunction of the vector
  // Laplacian.
  auto f = [&u](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 2 * pi * pi * u(x);
  };
  // Wrap both vectorfields into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  lf::mesh::utils::MeshFunctionGlobal mf_p(p);
  // ********** Part II: Loop over sequence of meshes **********
  // Generate a small unstructured triangular mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);

  // Number of levels
  const int L = multi_mesh.NumLevels();
  // Table of error norms
  std::vector<std::tuple<size_t, double, double>> L2errs;
  for (int level = 0; level < L; ++level) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // ********** Part III: Solving on a single mesh **********
    // Initialize dof handler for monolithic Whitney FEM
    lf::assemble::UniformFEDofHandler dofh(lev_mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 0}});
    const size_t N = dofh.NumDofs();
    LF_ASSERT_MSG(
        (N == lev_mesh_p->NumEntities(2) + lev_mesh_p->NumEntities(1)),
        "No dof mismatch");
    Eigen::VectorXd dofvec = solveMagStat2DBVP(dofh, mf_one, mf_f);
    // Build MeshFunction representing the 1-form component of the solution
    const MeshFunctionWF1 mf_sol_u(dofh, dofvec);
    const MeshFunctionWF0 mf_sol_p(dofh, dofvec);
    // Compute L2 norm of the error
    double L2err_u = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_u - mf_u), 2));
    double L2err_p = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_p - mf_p), 2));
    L2errs.push_back({N, L2err_u, L2err_p});
    std::cout << "L2 errors on level " << level << " : of u  = " << L2err_u
              << ", of p = " << L2err_p << std::endl;
    if (level > 0) {
      EXPECT_TRUE(get<1>(L2errs[level - 1]) > 1.9 * get<1>(L2errs[level]))
          << "Implementation does not reach expected convergence rate for u.";
    }
  }
}

}  // namespace MagStat2D::test
