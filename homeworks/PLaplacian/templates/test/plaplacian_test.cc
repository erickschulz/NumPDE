/**
 * @file plaplacian_test.cc
 * @brief NPDE homework plaplacian code
 * @author W. Tonnon
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../plaplacian.h"

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>

namespace plaplacian::test {

TEST(plaplacian, fixedPointSolvePLaplacian) {
  // Coefficient p for the p-Laplacian
  double p = 3;
  // Number of refinements in the mesh for estimating the order
  int reflevels = 4;

  // Manufactured solution of -div( ||grad u||^(p-2) grad u) = f with vanishing
  // Dirichlet BCs
  auto u = [](Eigen::Vector2d x) -> double {
    return sin(M_PI * x[0]) * sin(M_PI * x[1]);
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI),
                           M_PI * sin(x[0] * M_PI) * cos(x[1] * M_PI));
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return 2 * M_PI * M_PI * M_PI *
           (-8 * cos(x[0] * M_PI) * cos(x[0] * M_PI) * cos(x[1] * M_PI) *
                cos(x[1] * M_PI) +
            3 * cos(x[0] * M_PI) * cos(x[0] * M_PI) +
            3 * cos(x[1] * M_PI) * cos(x[1] * M_PI)) *
           sin(x[0] * M_PI) * sin(x[1] * M_PI) /
           sqrt(-cos(M_PI * (2 * x[0] - 2 * x[1])) -
                cos(M_PI * (2 * x[0] + 2 * x[1])) + 2);
  };
  // Wrap lambdas into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Triangular mesh hierarachy of unit square for testing
  // Adapted from Lehrfem++ \cppfile{homDir_linfe_demo.cc}
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << "\t Sequence of nested meshes used in test routine\n";
  multi_mesh.PrintInfo(std::cout);
  std::size_t L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeping error norms
  std::vector<std::tuple<int, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute finite-element solution of boundary value problem
    const Eigen::VectorXd sol_vec = fixedPointSolvePLaplacian(fe_space, f, p);

    // Compute error norms
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(mesh_p->NumEntities(2), L2err, H1serr);
  }
  // Output table of errors and estimated orders
  std::cout << "\t Table of error norms, p = " << p << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::left << std::setw(16)
            << "L2 error" << std::left << std::setw(16) << "L2 order"
            << std::left << std::setw(16) << "H1 error" << std::left
            << std::setw(16) << "H1 order" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int i = 1; i < reflevels + 1; ++i) {
    auto [N, l2err, h1serr] = errs.at(i);
    auto [Nprev, l2errprev, h1serrprev] = errs.at(i - 1);

    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::left << std::setw(16)
              << std::log(l2errprev / l2err) / std::log(2) << std::left
              << std::setw(16) << h1serr << std::left << std::setw(16)
              << std::log(h1serrprev / h1serr) / std::log(2) << std::endl;

    if (i >= 3) {
      EXPECT_GT(std::log(l2errprev / l2err) / std::log(2), 1.95);
      EXPECT_GT(std::log(h1serrprev / h1serr) / std::log(2), 0.95);
    }
  }
};

TEST(plaplacian, fixedPointSolvePLaplacian2) {
  // Coefficient p for the p-Laplacian
  double p = 1.5;
  // Number of refinements in the mesh for estimating the order
  int reflevels = 4;

  // Manufactured solution of -div( ||grad u||^(p-2) grad u) = f with vanishing
  // Dirichlet BCs
  auto u = [](Eigen::Vector2d x) -> double {
    return sin(M_PI * x[0]) * sin(M_PI * x[1]);
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI),
                           M_PI * sin(x[0] * M_PI) * cos(x[1] * M_PI));
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return -M_PI *
               (0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[1] * M_PI) *
                    sin(x[1] * M_PI) * cos(x[0] * M_PI) -
                0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * cos(x[0] * M_PI) *
                    cos(x[1] * M_PI) * cos(x[1] * M_PI)) *
               sin(x[1] * M_PI) * cos(x[0] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        1.25) -
           M_PI *
               (0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                    sin(x[1] * M_PI) * cos(x[1] * M_PI) -
                0.5 * M_PI * M_PI * M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI) *
                    cos(x[0] * M_PI) * cos(x[1] * M_PI)) *
               sin(x[0] * M_PI) * cos(x[1] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        1.25) +
           2 * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[1] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        0.25);
  };
  // Wrap lambdas into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Triangular mesh hierarachy of unit square for testing
  // Adapted from Lehrfem++ \cppfile{homDir_linfe_demo.cc}
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << "\t Sequence of nested meshes used in test routine\n";
  multi_mesh.PrintInfo(std::cout);
  std::size_t L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeping error norms
  std::vector<std::tuple<int, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute finite-element solution of boundary value problem
    const Eigen::VectorXd sol_vec = fixedPointSolvePLaplacian(fe_space, f, p);

    // Compute error norms
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(mesh_p->NumEntities(2), L2err, H1serr);
  }
  // Output table of errors and estimated orders
  std::cout << "\t Table of error norms, p = " << p << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::left << std::setw(16)
            << "L2 error" << std::left << std::setw(16) << "L2 order"
            << std::left << std::setw(16) << "H1 error" << std::left
            << std::setw(16) << "H1 order" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int i = 1; i < reflevels + 1; ++i) {
    auto [N, l2err, h1serr] = errs.at(i);
    auto [Nprev, l2errprev, h1serrprev] = errs.at(i - 1);

    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::setw(16)
              << std::log(l2errprev / l2err) / std::log(2) << std::setw(16)
              << h1serr << std::setw(16)
              << std::log(h1serrprev / h1serr) / std::log(2) << std::endl;

    if (i >= 3) {
      EXPECT_GT(std::log(l2errprev / l2err) / std::log(2), 1.9);
      EXPECT_GT(std::log(h1serrprev / h1serr) / std::log(2), 0.95);
    }
  }
};

TEST(plaplacian, newtonSolvePLaplacian) {
  // Coefficient p for the p-Laplacian
  double p = 3;
  // Number of refinements in the mesh for estimating the order
  int reflevels = 4;

  // Manufactured solution of -div( ||grad u||^(p-2) grad u) = f with vanishing
  // Dirichlet BCs
  auto u = [](Eigen::Vector2d x) -> double {
    return sin(M_PI * x[0]) * sin(M_PI * x[1]);
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI),
                           M_PI * sin(x[0] * M_PI) * cos(x[1] * M_PI));
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return 2 * M_PI * M_PI * M_PI *
           (-8 * cos(x[0] * M_PI) * cos(x[0] * M_PI) * cos(x[1] * M_PI) *
                cos(x[1] * M_PI) +
            3 * cos(x[0] * M_PI) * cos(x[0] * M_PI) +
            3 * cos(x[1] * M_PI) * cos(x[1] * M_PI)) *
           sin(x[0] * M_PI) * sin(x[1] * M_PI) /
           sqrt(-cos(M_PI * (2 * x[0] - 2 * x[1])) -
                cos(M_PI * (2 * x[0] + 2 * x[1])) + 2);
  };
  // Wrap lambdas into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Triangular mesh hierarachy of unit square for testing
  // Adapted from Lehrfem++ \cppfile{homDir_linfe_demo.cc}
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << "\t Sequence of nested meshes used in test routine\n";
  multi_mesh.PrintInfo(std::cout);
  std::size_t L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeM_PIng error norms
  std::vector<std::tuple<int, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute finite-element solution of boundary value problem
    const Eigen::VectorXd sol_vec = newtonSolvePLaplacian(fe_space, f, p);

    // Compute error norms
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(mesh_p->NumEntities(2), L2err, H1serr);
  }
  // Output table of errors and estimated orders
  std::cout << "\t Table of error norms, p = " << p << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::left << std::setw(16)
            << "L2 error" << std::left << std::setw(16) << "L2 order"
            << std::left << std::setw(16) << "H1 error" << std::left
            << std::setw(16) << "H1 order" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int i = 1; i < reflevels + 1; ++i) {
    auto [N, l2err, h1serr] = errs.at(i);
    auto [Nprev, l2errprev, h1serrprev] = errs.at(i - 1);

    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::setw(16)
              << std::log(l2errprev / l2err) / std::log(2) << std::setw(16)
              << h1serr << std::setw(16)
              << std::log(h1serrprev / h1serr) / std::log(2) << std::endl;

    if (i >= 3) {
      EXPECT_GT(std::log(l2errprev / l2err) / std::log(2), 1.95);
      EXPECT_GT(std::log(h1serrprev / h1serr) / std::log(2), 0.95);
    }
  }
};

TEST(plaplacian, newtonSolvePLaplacian2) {
  // Coefficient p for the p-Laplacian
  double p = 1.5;
  // Number of refinements in the mesh for estimating the order
  int reflevels = 4;

  // Manufactured solution of -div( ||grad u||^(p-2) grad u) = f with vanishing
  // Dirichlet BCs
  auto u = [](Eigen::Vector2d x) -> double {
    return sin(M_PI * x[0]) * sin(M_PI * x[1]);
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI),
                           M_PI * sin(x[0] * M_PI) * cos(x[1] * M_PI));
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return -M_PI *
               (0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[1] * M_PI) *
                    sin(x[1] * M_PI) * cos(x[0] * M_PI) -
                0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * cos(x[0] * M_PI) *
                    cos(x[1] * M_PI) * cos(x[1] * M_PI)) *
               sin(x[1] * M_PI) * cos(x[0] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        1.25) -
           M_PI *
               (0.5 * M_PI * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                    sin(x[1] * M_PI) * cos(x[1] * M_PI) -
                0.5 * M_PI * M_PI * M_PI * sin(x[1] * M_PI) * cos(x[0] * M_PI) *
                    cos(x[0] * M_PI) * cos(x[1] * M_PI)) *
               sin(x[0] * M_PI) * cos(x[1] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        1.25) +
           2 * M_PI * M_PI * sin(x[0] * M_PI) * sin(x[1] * M_PI) /
               std::pow(M_PI * M_PI * sin(x[0] * M_PI) * sin(x[0] * M_PI) *
                                cos(x[1] * M_PI) * cos(x[1] * M_PI) +
                            M_PI * M_PI * sin(x[1] * M_PI) * sin(x[1] * M_PI) *
                                cos(x[0] * M_PI) * cos(x[0] * M_PI),
                        0.25);
  };
  // Wrap lambdas into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Triangular mesh hierarachy of unit square for testing
  // Adapted from Lehrfem++ \cppfile{homDir_linfe_demo.cc}
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << "\t Sequence of nested meshes used in test routine\n";
  multi_mesh.PrintInfo(std::cout);
  std::size_t L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeping error norms
  std::vector<std::tuple<int, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute finite-element solution of boundary value problem
    const Eigen::VectorXd sol_vec = newtonSolvePLaplacian(fe_space, f, p);

    // Compute error norms
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(mesh_p->NumEntities(2), L2err, H1serr);
  }
  // Output table of errors and estimated orders
  std::cout << "\t Table of error norms, p = " << p << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::left << std::setw(16)
            << "L2 error" << std::left << std::setw(16) << "L2 order"
            << std::left << std::setw(16) << "H1 error" << std::left
            << std::setw(16) << "H1 order" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int i = 1; i < reflevels + 1; ++i) {
    auto [N, l2err, h1serr] = errs.at(i);
    auto [Nprev, l2errprev, h1serrprev] = errs.at(i - 1);

    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::setw(16)
              << std::log(l2errprev / l2err) / std::log(2) << std::setw(16)
              << h1serr << std::setw(16)
              << std::log(h1serrprev / h1serr) / std::log(2) << std::endl;

    if (i >= 3) {
      EXPECT_GT(std::log(l2errprev / l2err) / std::log(2), 1.9);
      EXPECT_GT(std::log(h1serrprev / h1serr) / std::log(2), 0.95);
    }
  }
};

}  // namespace plaplacian::test
