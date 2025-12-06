/**
 * @file LeastSquaresAdvection_test.cc
 * @brief NPDE homework LeastSquaresAdvection code
 * @author Ralf Hiptmair
 * @date July 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../leastsquaresadvection.h"

#include <Eigen/src/Core/Matrix.h>
#include <gtest/gtest.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/fe/fe_tools.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_global.h>

#include <Eigen/Core>

/* Tests in the google testing framework

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

namespace LeastSquaresAdvection::test {


TEST(LeastSquaresAdvection, SUBPROBLEM_F) {
  // Build "mesh" comprising a single triangle
  // Short name for 2d coordinate vectors
  using coord_t = Eigen::Vector2d;
  // Corner coordinates for a triangle
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;
  // Ftehc a mesh facvtory
  const std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Add points
  mesh_factory_ptr->AddPoint(coord_t({0, 0}));    // point 0
  mesh_factory_ptr->AddPoint(coord_t({1, 0.2}));  // point 1
  mesh_factory_ptr->AddPoint(coord_t({0.3, 1}));  // point 2
  // Define triangular cells
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({0, 1, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Ready to build the mesh data structure
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();
  // Pointer to single cell of the mesh
  const lf::mesh::Entity* cell_p = mesh_p->EntityByIndex(0, 0);
  LF_ASSERT_MSG(cell_p != nullptr, "Invalid cell!");
  // Mesh function for reaction coefficient
  lf::mesh::utils::MeshFunctionGlobal mf_kappa{
      [](Eigen::Vector2d x) -> double { return x[0] * x[1]; }};
  // Velocity vector
  Eigen::Vector2d velocity(2.0, 1.0);

  constexpr bool generate_matrices = false;

  {
    // Linear Lagrangian FE
    // Construct finite element space
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Initialize element matrix provider object
    LeastSquaresAdvection::LSQAdvectionMatrixProvider emp(fe_space, velocity,
                                                          mf_kappa);
    // Obtain element matrix
    auto elmat = emp.Eval(*cell_p);
    if (generate_matrices) {
      std::cout << "Element matrix for p = 1:\n" << elmat << std::endl;
    } else {
      Eigen::MatrixXd elmat_ref(3, 3);
      elmat_ref << 2.73203, -2.12232, -0.804009, -2.12232, 1.6503, 0.624527,
          -0.804009, 0.624527, 0.239119;

      EXPECT_NEAR((elmat - elmat_ref).norm(), 0.0, 1.0E-3);
    }
  }
  {
    // Quadaratic Lagrangian FE
    // Construct finite element space
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
    // Initialize element matrix provider object
    LeastSquaresAdvection::LSQAdvectionMatrixProvider emp(fe_space, velocity,
                                                          mf_kappa);
    // Obtain element matrix
    auto elmat = emp.Eval(*cell_p);
    if (generate_matrices) {
      std::cout << "Element matrix for p = 2:\n" << elmat << std::endl;
    } else {
      Eigen::MatrixXd elmat_ref(6, 6);
      elmat_ref << 2.80047, 0.69486, 0.25289, -2.76202, 0.0152392, -0.963411,
          0.69486, 1.58402, -0.200868, -2.78049, 0.766799, 0.0121759, 0.25289,
          -0.200868, 0.213577, 0.00998669, 0.775315, -1.0147, -2.76202,
          -2.78049, 0.00998669, 5.95794, -2.0079, 1.40248, 0.0152392, 0.766799,
          0.775315, -2.0079, 6.38509, -5.60093, -0.963411, 0.0121759, -1.0147,
          1.40248, -5.60093, 5.87888;
      EXPECT_NEAR((elmat - elmat_ref).norm(), 0.0, 1.0E-3);
    }
  }
  {
    // Cubic Lagrangian FE
    // Construct finite element space
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO3<double>>(mesh_p);
    // Initialize element matrix provider object
    LeastSquaresAdvection::LSQAdvectionMatrixProvider emp(fe_space, velocity,
                                                          mf_kappa);
    // Obtain element matrix
    auto elmat = emp.Eval(*cell_p);
    if (generate_matrices) {
      std::cout << "Element matrix for p = 3:\n" << elmat << std::endl;
    }
  }
}

TEST(LeastSquaresAdvection, SUBPROBLEM_G) {
  // Simple mesh of the unit square
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Construct finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Mesh function for reaction coefficient = 0
  lf::mesh::utils::MeshFunctionGlobal mf_kappa{
      [](Eigen::Vector2d /*x*/) -> double { return 0.0; }};
  // Velocity vector
  Eigen::Vector2d velocity(2.0, 1.0);
  // Mesh function for boundary data = 1
  lf::mesh::utils::MeshFunctionGlobal mf_g{
      [](Eigen::Vector2d x) -> double { return sin(x(0)); }};
  Eigen::VectorXd mu_vec = LeastSquaresAdvection::solveAdvectionDirichletBVP(
      fe_space, velocity, mf_kappa, mf_g);
  ASSERT_EQ(mu_vec.size(), 40);
  Eigen::Vector<double, 40> mu_vec_exact(mu_vec.size());
  mu_vec_exact << 0, 0.479426, 0.841471, -0.0220163, 0.0756852, 0.341568, 0,
      -0.0263345, 0.00251743, -0.0263345, 0, -0.00121418, 0.0032156, 0.247404,
      -0.0132169, 0, 0.681639, 0.127439, 0.257852, 0.471159, 0.59942,
      0.00689166, -0.00226944, -0.0241754, -0.00381867, 0.202784, 0.00703468,
      0.0873403, 0.0873403, 0.00125872, 0, -0.00536789, -0.0263345, 0.00173787,
      -0.00649425, 0.000187518, 0.00234808, -0.00649425, -0.00106241, 0.0038499;

  // Exact solution should = 1
  EXPECT_NEAR((mu_vec - mu_vec_exact).norm(), 0.0, 1E-6);
}

}  // namespace LeastSquaresAdvection::test

/*
  Element matrix for p = 1:
  2.73203  -2.12232 -0.804009
 -2.12232    1.6503  0.624527
-0.804009  0.624527  0.239119
Element matrix for p = 2:
   2.80047    0.69486    0.25289   -2.76202  0.0152392  -0.963411
   0.69486    1.58402  -0.200868   -2.78049   0.766799  0.0121759
   0.25289  -0.200868   0.213577 0.00998669   0.775315    -1.0147
  -2.76202   -2.78049 0.00998669    5.95794    -2.0079    1.40248
 0.0152392   0.766799   0.775315    -2.0079    6.38509   -5.60093
 -0.963411  0.0121759    -1.0147    1.40248   -5.60093    5.87888
Element matrix for p = 3:
    2.39089   -0.365629   -0.131899    -2.75344     1.45841      -0.214
-0.215167     0.64762   -0.833809 -0.00205298 -0.365629     1.33484
0.101846     1.36849    -2.85902      0.9208    -0.22338   -0.119329   -0.116109
-0.00188096 -0.131899    0.101846    0.176283   -0.015083   -0.016896 -0.328938
0.830291    -1.06342    0.469021 -0.00177148 -2.75344     1.36849
-0.015083     7.64158    -5.73431     0.50722    0.501328 -0.361207      1.8208
-2.97675 1.45841    -2.85902   -0.016896    -5.73431     7.60743    -2.49959
0.505934   -0.357793   -0.360922     2.18207 -0.214      0.9208   -0.328938
0.50722    -2.49959     7.90966   -0.440966     1.41882      1.4081    -8.44098
  -0.215167    -0.22338    0.830291    0.501328    0.505934
-0.440966     7.93208    -7.04949      1.4191    -2.99804 0.64762   -0.119329
-1.06342   -0.361207   -0.357793     1.41882    -7.04949     7.51361 -2.99468
2.17273 -0.833809   -0.116109    0.469021      1.8208
-0.360922      1.4081      1.4191    -2.99468     7.61446    -8.42758
-0.00205298 -0.00188096 -0.00177148    -2.97675     2.18207    -8.44098 -2.99804
2.17273    -8.42758     18.2412
*/
