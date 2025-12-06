/**
 * @file IRKDegenerateEvl_test.cc
 * @brief NPDE homework IRKDegenerateEvl code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../irkdegenerateevl.h"

#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <gtest/gtest.h>
#include <lf/assemble/coomatrix.h>
#include <lf/fe/fe_tools.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

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

namespace IRKDegenerateEvl::test {

// TEST(IRKDegenerateEvl, PleaseNameTest) { EXPECT_TRUE(true); }

TEST(IRKDegenerateEvl, TestM) {
  // Simple mesh of the unit square
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Construct finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Compute basis expansion coefficients of a quadratic function
  lf::mesh::utils::MeshFunctionGlobal mf_quad(
      [](Eigen::Vector2d x) -> double { return (x[0] * x[1] + 1.0); });
  Eigen::VectorXd mu_vec = lf::fe::NodalProjection(*fe_space, mf_quad);
  // Assemble boundary mass matrix
  lf::assemble::COOMatrix<double> M{IRKDegenerateEvl::buildM(fe_space)};
  // Compute squared L2 norm of the function on the boundary
  const double res = mu_vec.dot(M.MatVecMult(1.0, mu_vec));

  std::cout << res << std::endl;
  EXPECT_NEAR(res, 20.0 / 3.0, 1.0E-8);
}


TEST(IRKDegenerateEvl, SUBPROBLEM_E) {
  // Simple mesh of the unit square
  // Short name for 2d coordinate vectors
  using coord_t = Eigen::Vector2d;
  // Corner coordinates for a triangle
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;
  // Ftehc a mesh facvtory
  const std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Add points
  mesh_factory_ptr->AddPoint(coord_t({0, 0}));        // point 0
  mesh_factory_ptr->AddPoint(coord_t({1, 0}));        // point 1
  mesh_factory_ptr->AddPoint(coord_t({1, 1}));        // point 2
  mesh_factory_ptr->AddPoint(coord_t({0, 1}));        // point 3
  mesh_factory_ptr->AddPoint(coord_t({0.31, 0.27}));  // point 4

  // Define triangular cells
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({0, 1, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({1, 2, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({2, 3, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({3, 0, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Ready to build the mesh data structure
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // Construct finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Assemble boundary mass matrix M
  const lf::assemble::COOMatrix<double> M{IRKDegenerateEvl::buildM(fe_space)};
  // Convert into a dense matrix
  const Eigen::MatrixXd M_dense{M.makeDense()};

  constexpr bool generate_matrix = false;
  if (generate_matrix) {
    // define the format you want, you only need one instance of this...
    // see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                           Eigen::DontAlignCols, ", ", ",\n");
    std::cout << "M = \n " << M_dense.format(CSVFormat) << std::endl;
  } else {
    Eigen::MatrixXd M_ref(M.rows(), M.cols());
    M_ref << 0.266667, -0.0333333, 0, -0.0333333, 0, 0.0666667, 0.0666667, 0, 0,
        0, 0, 0, 0, -0.0333333, 0.266667, -0.0333333, 0, 0, 0.0666667, 0, 0,
        0.0666667, 0, 0, 0, 0, 0, -0.0333333, 0.266667, -0.0333333, 0, 0, 0, 0,
        0.0666667, 0, 0.0666667, 0, 0, -0.0333333, 0, -0.0333333, 0.266667, 0,
        0, 0.0666667, 0, 0, 0, 0.0666667, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0.0666667, 0.0666667, 0, 0, 0, 0.533333, 0, 0, 0, 0, 0, 0, 0,
        0.0666667, 0, 0, 0.0666667, 0, 0, 0.533333, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0666667, 0.0666667, 0, 0, 0, 0, 0,
        0.533333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0666667, 0.0666667, 0, 0, 0, 0, 0, 0, 0.533333, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    EXPECT_NEAR((M_dense - M_ref).norm(), 0.0, 1.0E-3);
  }
}

template <typename MATRIX>
void testTimestepping(const MATRIX& Ark, unsigned int n_ref = 7,
                      unsigned int min_steps = 10) {
  // Simple mesh of the unit square
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Construct finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Build MOL ODE matrices
  lf::assemble::COOMatrix<double> M{IRKDegenerateEvl::buildM(fe_space)};
  lf::assemble::COOMatrix<double> A{IRKDegenerateEvl::buildA(fe_space)};
  // Initial vector
  const Eigen::VectorXd mu0 = Eigen::VectorXd::LinSpaced(A.cols(), 0.0, 1.0);
  // Run timestepping with different number of timesteps
  std::vector<Eigen::VectorXd> mu_fin(n_ref + 1, Eigen::VectorXd(A.cols()));
  unsigned int n_ts = min_steps;
  const double mu0_avg =
      mu0.dot(M.MatVecMult(1.0, Eigen::VectorXd::Constant(A.cols(), 1.0)));
  for (int l = 0; l <= n_ref; ++l, n_ts *= 2) {
    mu_fin[l] =
        IRKDegenerateEvl::timesteppingIRKMOLODE(M, A, Ark, mu0, n_ts, 1.0);
    const double mu_fin_avg = mu_fin[l].dot(
        M.MatVecMult(1.0, Eigen::VectorXd::Constant(A.cols(), 1.0)));
    EXPECT_NEAR(mu_fin_avg, mu0_avg, 1.0E-8);
  }
  n_ts = min_steps;
  std::cout << Ark.cols() << "-stage IRK SSM with Butcher matrix\n"
            << Ark << std::endl;

  double err_old;
  for (int l = 0; l < n_ref; ++l, n_ts *= 2) {
    double err = (mu_fin[l] - mu_fin[n_ref]).norm();
    std::cout << n_ts << " steps: error = " << err;
    if (l > 0) {
      double ratio = err_old / err;
      std::cout << ", ratio = " << ratio;
      if (l > 4) {
        if (Ark.cols() == 1)
          ASSERT_GE(ratio, 1.99) << "Convergence rate for 1-stage SSM too low";
        if (Ark.cols() == 2)
          ASSERT_GE(ratio, 7.98) << "Convergence rate for 2-stage SSM too low";
      }
    }
    std::cout << std::endl;
    err_old = err;
  }
}

TEST(IRKDegenerateEvl, Ieul) {
  // Butcher "matrix"
  const Eigen::Matrix<double, 1, 1> Ark =
      (Eigen::Matrix<double, 1, 1>() << 1).finished();
  testTimestepping(Ark);
}

TEST(IRKDegenerateEvl, Radau) {
  // Butcher "matrix"
  const Eigen::Matrix<double, 2, 2> Ark =
      (Eigen::Matrix<double, 2, 2>() << 5.0 / 12.0, -1.0 / 12.0, 3.0 / 4.0,
       1.0 / 4.0)
          .finished();
  testTimestepping(Ark);
}

TEST(IRKDegenerateEvl, Tabulate) {
  // Simple mesh of the unit square
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Construct finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  double T = 1.0;
  unsigned int no_ts = 10;
  auto it_norms = IRKDegenerateEvl::tabulateSolNorms(
      fe_space,
      [](Eigen::Vector2d x) -> double {
        return std::cos(0.5 * M_PI * x[0]) * std::cos(0.5 * M_PI * x[1]);
      },
      no_ts, T);
  size_t m = it_norms.size();
  ASSERT_GE(m, 2);
  EXPECT_GE(it_norms[0].first, it_norms[m - 1].first)
      << "L2-nNorm should decrease!";
  EXPECT_GE(it_norms[0].second, it_norms[m - 1].second)
      << "H1-norm should decrease!";
}


TEST(IRKDegenerateEvl, SUBPROBLEM_F) {
  // Define the Butcher Tableau matrix (Ark)
  Eigen::MatrixXd Ark(5, 5);
  Ark << 0., 0., 0., 0., 0., 1. / 3., 1. / 3., 0., 0., 0., 2. / 3., -1. / 3.,
      1., 0., 0., 1., 1., -1., 1., 0., 0., 1. / 8., 3. / 8., 3. / 8., 1. / 8.;

  // Define the initial vector (mu0)
  Eigen::VectorXd mu0(5);
  mu0 << 1., 1., 1., 1., 1.;

  // Create and initialize the matrices M and A
  lf::assemble::COOMatrix<double> M(5, 5);
  lf::assemble::COOMatrix<double> A(5, 5);
  for (int i = 0; i < 5; ++i) M.AddToEntry(i, i, 1);
  for (int i = 0; i < 5; ++i) A.AddToEntry(i, i, i + 1);

  // Time-stepping parameters
  unsigned int no_ts = 250;  // Number of timesteps
  double T = 1.;             // Total time

  // Call the IRK MOLODE time-stepping function
  auto out = IRKDegenerateEvl::timesteppingIRKMOLODE(M, A, Ark, mu0, no_ts, T);

  // Validate the results using EXPECT_NEAR
  for (int i = 0; i < 5; ++i) {
    EXPECT_NEAR(out(i), exp(-i - 1), 1e-2);
  }
}

TEST(IRKDegenerateEvl, SUBPROBLEM_H) {
  // Create a mock finite element space
  // Simple mesh of the unit square
  // Short name for 2d coordinate vectors
  using coord_t = Eigen::Vector2d;
  // Corner coordinates for a triangle
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;
  // Ftehc a mesh facvtory
  const std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Add points
  mesh_factory_ptr->AddPoint(coord_t({0, 0}));        // point 0
  mesh_factory_ptr->AddPoint(coord_t({1, 0}));        // point 1
  mesh_factory_ptr->AddPoint(coord_t({1, 1}));        // point 2
  mesh_factory_ptr->AddPoint(coord_t({0, 1}));        // point 3
  mesh_factory_ptr->AddPoint(coord_t({0.31, 0.27}));  // point 4

  // Define triangular cells
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({0, 1, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({1, 2, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({2, 3, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({3, 0, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Ready to build the mesh data structure
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  auto fes = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Define a simple initial condition: u0(x, y) = x + y
  auto u0 = [](const Eigen::Vector2d& x) -> double { return x[0] + x[1]; };

  // Define time-stepping parameters
  unsigned int no_ts = 30;  // Number of time steps
  double T = 10.0;          // Final time

  // Run the function
  auto norms = IRKDegenerateEvl::tabulateSolNorms(fes, u0, no_ts, T);

  // Verify the computed norms (basic checks)
  ASSERT_EQ(norms.size(),
            no_ts + 1);  // There should be one norm pair per time step

  // Optionally print results for debugging purposes
  std::cout << "Time step norms (L2 and H1):\n";
  for (size_t i = 0; i < norms.size(); ++i) {
    std::cout << "Step " << i << ": L2 = " << norms[i].first
              << ", H1 = " << norms[i].second << "\n";
  }

  // Specific norm checks (this requires understanding of the problem dynamics)
  EXPECT_NEAR(2., norms.back().first, 1e-5);
  EXPECT_NEAR(0., norms.back().second, 1e-5);
}

}  // namespace IRKDegenerateEvl::test
