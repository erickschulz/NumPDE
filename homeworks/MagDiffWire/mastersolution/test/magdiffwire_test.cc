/**
 * @file MagDiffWire_test.cc
 * @brief NPDE homework MagDiffWire code
 * @author Ralf Hiptmair
 * @date JUne 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../magdiffwire.h"

#include <gtest/gtest.h>
#include <lf/base/lf_assert.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>
#include <numbers>

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

namespace MagDiffWire::test {

TEST(MagDiffWire, MatTest) {
  using namespace std::numbers;
  using std::cos;
  using std::sin;

  // Function to be interpolated
  auto u = [](Eigen::VectorXd x) -> double {
    return sin(pi * x[0]) * sin(pi * x[1]);
  };
  // Wrap into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  // Coefficients
  const double sigma_c = 1.0;
  const double mu_c = 1.0;

  // Compute integrals on a sequence of regularly refined meshes
  unsigned int refsteps = 4;
  // Generate a small unstructured triangular mesh of the unit square
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Generate sequence of meshes by uniform regular refinement
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  // std::cout << "\t Sequence of nested meshes created\n";
  // multi_mesh.PrintInfo(std::cout);
  // Number of levels
  const int L = multi_mesh.NumLevels();
  // Integral values
  double u_int;
  double u_sq_int;
  double u_grad_sq_int;
  double u_int2;
  for (int level = 0; level < L; ++level) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Initialize Lagrangian finite element space
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO2<double>> fes_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);
    // Number of FE degrees of freedom
    const size_t N = fes_p->LocGlobMap().NumDofs();
    // Setting: constant coefficients = 1, all cells flagged as located in the
    // conducting domain
    const lf::mesh::utils::CodimMeshDataSet<bool> all_cells(lev_mesh_p, 0,
                                                            true);
    // Obtain extended Galerkin matrices as lf::asseemble::COOMatrix
    auto [Mt, At] = buildExtMOLMatrices(fes_p, all_cells, sigma_c, mu_c);
    LF_ASSERT_MSG(At.rows() == N + 1, "A wrong row number");
    LF_ASSERT_MSG(At.cols() == N + 1, "A wrong column number");
    LF_ASSERT_MSG(Mt.rows() == N + 1, "M wrong row number");
    LF_ASSERT_MSG(Mt.cols() == N + 1, "M wrong col number");
    // Computer coefficient vector for nodal interpolant
    Eigen::VectorXd mu_vec(N + 1);
    mu_vec.head(N) = lf::fe::NodalProjection(*fes_p, mf_u);
    mu_vec[N] = 0.0;
    // Compute various norms/integrals
    Eigen::VectorXd vec_one = Eigen::VectorXd::Constant(N + 1, 1.0);
    u_int = vec_one.dot(Mt.MatVecMult(1.0, mu_vec));
    u_sq_int = mu_vec.dot(Mt.MatVecMult(1.0, mu_vec));
    u_grad_sq_int = mu_vec.dot(At.MatVecMult(1.0, mu_vec));
    u_int2 = At.MatVecMult(1.0, mu_vec)[N];
    /*
    std::cout << "level " << level << ", N = " << N << ", u_int = " << u_int
              << ", u_sq_int = " << u_sq_int
              << ", u_grad_sq_int = " << u_grad_sq_int
              << ", u_int2 = " << u_int2 << std::endl;
    */
  }
  /*
  std::cout << "\nExact values\n"
            << ", u_int = " << 4.0 / (pi * pi) << ", u_sq_int = " << 0.25
            << ", u_grad_sq_int = " << (pi * pi / 2.0)
            << ", u_int2 = " << 4 / (pi * pi) << std::endl;
  */

  EXPECT_NEAR(u_int, 4.0 / (pi * pi), 0.01);
  EXPECT_NEAR(u_sq_int, 0.25, 0.01);
  EXPECT_NEAR(u_grad_sq_int, (pi * pi / 2.0), 0.01);
  EXPECT_NEAR(u_int2, 4 / (pi * pi), 0.01);
}

TEST(MagDiffWire, sdirkTest) {
  int refsteps = 3;
  using namespace std::numbers;
  using std::cos;
  using std::exp;
  using std::sin;
  // Manufactured solution on the unit square
  auto u = [](Eigen::VectorXd x, double t) -> double {
    return sin(pi * x[0]) * sin(pi * x[1]) * exp(-2 * pi * pi * t);
  };
  auto u0 = [&u](Eigen::VectorXd x) -> double { return u(x, 0.0); };
  auto grad_u = [](Eigen::VectorXd x, double t) -> Eigen::Vector2d {
    return pi * exp(-2 * pi * pi * t) *
           Eigen::Vector2d(cos(pi * x[0]) * sin(pi * x[1]),
                           sin(pi * x[0]) * cos(pi * x[1]));
  };
  lf::mesh::utils::MeshFunctionGlobal mf_u0(u0);
  const double T_final = 0.5;  // Final time T
                               // Solution at final time
  auto u_fin = [&u, T_final](Eigen::Vector2d x) -> double {
    return u(x, T_final);
  };
  auto grad_u_fin = [&grad_u, T_final](Eigen::VectorXd x) -> Eigen::Vector2d {
    return grad_u(x, T_final);
  };
  // Wrap into MeshFunctions
  lf::mesh::utils::MeshFunctionGlobal mf_u_fin(u_fin);
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u_fin(grad_u_fin);

  // Exciting current
  auto I_source = [](double t) -> double {
    return 4.0 / (pi * pi) * exp(-2 * pi * pi * t);
  };

  // Generate a small unstructured triangular mesh of the unit square
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Generate sequence of meshes by uniform regular refinement
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);

  // Number of levels
  const int L = multi_mesh.NumLevels();
  // Loop over the levels
  const unsigned int M0 = 20;  // No. of timesteps on coarsest mesh
  unsigned int M = M0;         // Number of timesteps
  // Table of various error norms
  std::vector<std::tuple<size_t, double, double>> errs;
  for (int level = 0; level < L; ++level, M *= 2) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Initialize Lagrangian finite element space
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO2<double>> fes_o2_ptr =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);
    // Number of FE degrees of freedom
    const size_t N = fes_o2_ptr->LocGlobMap().NumDofs();
    // Setting: constant coefficients = 1, all cells flagged as located in the
    // conducting domain
    const lf::mesh::utils::CodimMeshDataSet<bool> all_cells(lev_mesh_p, 0,
                                                            true);
    // Set initial data: nodal interpolant of u(x,0)
    Eigen::VectorXd mu_vec(N + 1);
    mu_vec.head(N) = lf::fe::NodalProjection(*fes_o2_ptr, mf_u0);
    mu_vec[N] = 0.0;
    // Main timestepping loop
    sdirkMagDiffWire(fes_o2_ptr, all_cells, 1.0, 1.0, I_source, M, T_final,
                     mu_vec);
    // Compute errors at final time
    const lf::fe::MeshFunctionFE mf_uh_fin(fes_o2_ptr, mu_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_uh_fin(fes_o2_ptr, mu_vec);
    const double L2_err = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_uh_fin - mf_u_fin), 4));
    const double H1s_err = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p,
        lf::mesh::utils::squaredNorm(mf_grad_uh_fin - mf_grad_u_fin), 4));
    errs.emplace_back(N, L2_err, H1s_err);
  }
  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
  std::cout.precision(3);
  std::cout << std::endl
            << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 err" << std::setw(16) << "H1s err" << std::setw(16)
            << "L2 ratio" << '\n';
  std::cout << "---------------------------------------------" << '\n';
  double L2err_prev = 100.;
  for (const auto& err : errs) {
    auto [N, L2err, H1err] = err;
    double ratio = L2err_prev / L2err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << L2err << std::setw(16) << H1err << std::setw(16) << ratio
             << '\n';
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << L2err << std::setw(16) << H1err << std::setw(16) << ratio
              << '\n';
    L2err_prev = L2err;
    EXPECT_TRUE(ratio > 1.9);
  }
}

TEST(MagDiffWire, DAETest) {
  // Testing SDIRK-2 for an index-1 DAW
  const Eigen::Matrix2d M = (Eigen::Matrix2d() << 1, 0, 0, 0).finished();
  const Eigen::Matrix2d A = (Eigen::Matrix2d() << 2, -1, -1, 2).finished();
  auto rho = [](double t) -> Eigen::Vector2d {
    return Eigen::Vector2d(0.0, std::sin(t));
  };
  // Final time
  const double T_final = 10;
  const double zeta = 1.0 - 0.5 * std::sqrt(2.0);  // SDIRK parameter
  // Smallest number of timesteps
  unsigned int N = 20;
  // Number of refinements
  const unsigned int refsteps = 10 + 2;
  std::vector<std::pair<double, Eigen::Vector2d>> mu_vecs_fin{};
  for (unsigned int l = 0; l < refsteps; ++l, N *= 2) {
    // Initial value
    Eigen::Vector2d mu_vec(1.0, 0.0);
    // Timestep size
    const double tau = T_final / N;
    // Matrix "to be inverted"
    Eigen::MatrixXd D = M + tau * zeta * A;

    // Starting time
    double t = 0.0;
    for (int k = 1; k <= N; ++k, t += tau) {
      const Eigen::Vector2d Amu = A * mu_vec;
      Eigen::Vector2d rho_vec = rho(t + zeta * tau);
      const Eigen::Vector2d kappa1 = D.lu().solve(rho_vec - Amu);
      rho_vec = rho(t + tau);
      const Eigen::Vector2d kappa2 =
          D.lu().solve(rho_vec - tau * (1 - zeta) * A * kappa1 - Amu);
      mu_vec += tau * (1 - zeta) * kappa1 + tau * zeta * kappa2;
    }
    mu_vecs_fin.push_back({tau, mu_vec});
  }
  double err, err_old = 0.0;
  for (unsigned int l = 4; l < refsteps - 2; l++) {
    std::cout << "tau = " << mu_vecs_fin[l].first << ", err = "
              << (err =
                      (mu_vecs_fin[l].second - mu_vecs_fin[refsteps - 1].second)
                          .norm());
    if (err_old > 0) {
      std::cout << ", ratio = " << err_old / err;
      std::cout << std::endl;
      EXPECT_TRUE(err_old > 4. * err);
    }
    err_old = err;
  }
}

}  // namespace MagDiffWire::test
