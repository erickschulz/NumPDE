/**
 * @file minimalgraphsurface_test.cc
 * @brief NPDE homework "MinimalGraphSurface" code
 * @author Wouter Tonnon
 * @date 12.04.2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../minimalgraphsurface.h"

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/fe/fe_tools.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_traits.h>
#include <lf/uscalfe/fe_space_lagrange_o1.h>

#include <Eigen/Core>

namespace MinimalGraphSurface::test {

TEST(MinimalGraphSurface, computeGraphArea) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is simple, so we can compute the graph area
  // analytically
  auto mf = lf::mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2. * x[1]; });

  // Project onto finite-element space
  auto dof_vec = lf::fe::NodalProjection(*fe_space, mf);
  auto mf_fe = lf::fe::MeshFunctionFE(fe_space, dof_vec);

  // Compute the surface are using computeGraphArea
  auto integral_val = MinimalGraphSurface::computeGraphArea(fe_space, dof_vec);

  // Test if the value is correct
  double tol = 1e-2;
  EXPECT_NEAR(integral_val, 9. * std::sqrt(6.), tol);
}

TEST(MinimalGraphSurface, CoeffScalarc) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is simple, so we can compute the coefficient c
  // analytically
  auto mf = lf::mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2. * x[1]; });

  // Project onto finite-element space
  auto dof_vec = lf::fe::NodalProjection(*fe_space, mf);
  auto mf_fe = lf::fe::MeshFunctionFE(fe_space, dof_vec);

  // Compute the coefficient c using CoeffScalarc
  MinimalGraphSurface::CoeffScalarc c(fe_space, dof_vec);

  // Integrate the coefficient over the entire domain
  double integral_val = lf::fe::IntegrateMeshFunction(*mesh_p, c, 2);

  // Test if the value is correct
  double tol = 1e-2;
  EXPECT_NEAR(integral_val, -9. / sqrt(6), tol);

  // Test if CoeffScalarc satisfies the conditions for being a MeshFunction.
  EXPECT_TRUE(lf::mesh::utils::MeshFunction<MinimalGraphSurface::CoeffScalarc>);
}

TEST(MinimalGraphSurface, CoeffTensorA) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is simple, so we can compute the coefficient A
  // analytically
  auto mf = lf::mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2. * x[1]; });

  // Project onto finite-element space
  auto dof_vec = lf::fe::NodalProjection(*fe_space, mf);
  auto mf_fe = lf::fe::MeshFunctionFE(fe_space, dof_vec);

  // Compute the coefficient A using CoeffTensorA
  MinimalGraphSurface::CoeffTensorA Af(fe_space, dof_vec);

  // For the chosen function, the solution should be constant, so it does not
  // matter where we evaluate the coefficient. We choose to evaluate the
  // coefficient in local coordinates (0.25,0.25) corresponding to entity 0
  // (with codim 0).
  auto e = mesh_p->EntityByIndex(0, 0);
  std::vector<Eigen::Matrix<double, 2, 2>> A_mat =
      Af(*e, Eigen::Matrix<double, 2, 1>({.25, .25}));

  // The exact solution is 1/(6*sqrt(6))*[5.5, 5; 5, 4]

  double tol = 1e-8;
  EXPECT_NEAR(A_mat[0](0, 0), 1. / sqrt(6.) * (1. - 1. / 6. * 1.), tol);
  EXPECT_NEAR(A_mat[0](1, 0), 1. / sqrt(6.) * (-1. / 6. * 2.), tol);
  EXPECT_NEAR(A_mat[0](0, 1), 1. / sqrt(6.) * (-1. / 6. * 2.), tol);
  EXPECT_NEAR(A_mat[0](1, 1), 1. / sqrt(6.) * (1. - 1. / 6. * 4.), tol);

  // Test if CoeffTensorA satisfies the conditions for being a MeshFunction.
  EXPECT_TRUE(lf::mesh::utils::MeshFunction<MinimalGraphSurface::CoeffTensorA>);
}

TEST(MinimalGraphSurface, computeNewtonCorrection) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is linear, so we know that it is a minimal surface
  auto mf = lf::mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return 2 * x[0] - 3 * x[1]; });

  // Project onto finite-element space
  auto dof_vec = lf::fe::NodalProjection(*fe_space, mf);
  auto mf_fe = lf::fe::MeshFunctionFE(fe_space, dof_vec);

  // Compute the Newton correction using computerNewtonCorrection
  Eigen::VectorXd mu_vec =
      MinimalGraphSurface::computeNewtonCorrection(fe_space, dof_vec);

  // Since mf already is a minimal surface, the Newton correction should be zero
  double tol = 1e-8;
  EXPECT_NEAR(mu_vec.norm(), 0., tol);
}

TEST(MinimalGraphSurface, computeNewtonCorrection2) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_0_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // The generated mesh is too course to resolve the perturbation we will add
  // later. Therefor, we refine the mesh.
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_0_p, 4);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is linear, so we know that it is a minimal surface
  auto min_surface = lf::mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return 2 * x[0] - 3 * x[1]; });

  // We take a surface that is slightly perturbed compared to the minimal
  // surface
  auto perturbed_min_surface =
      lf::mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return 2 * x[0] - 3 * x[1] +
               .1 * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
      });

  // Project the perturbed minimal surface onto the finite-element space
  auto dof_vec = lf::fe::NodalProjection(*fe_space, perturbed_min_surface);

  // We now start with the perturbed minimal surface as initial condition and
  // check if it converges to the minimal surface when the Newton corrections
  // are used as corrections.
  int i = 0;
  Eigen::VectorXd mu_vec;
  do {
    mu_vec = MinimalGraphSurface::computeNewtonCorrection(fe_space, dof_vec);
    dof_vec += mu_vec;
  } while (mu_vec.norm() > 1e-10 && ++i < 1000);

  // We generate global mesh functions to compute the error
  lf::mesh::utils::MeshFunctionGlobal exact_mf{min_surface};
  auto mu_mf = lf::fe::MeshFunctionFE(fe_space, dof_vec);

  // The perturbed minimal surface should have converged to the minimal surface
  double L2err = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(exact_mf - mu_mf), 2));
  double tol = 1e-8;
  EXPECT_NEAR(L2err, 0., tol);
}

TEST(MinimalGraphSurface, graphMinimalSurface1) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_0_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // The generated mesh is too course . Therefore, we refine the mesh.
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_0_p, 4);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Our prototype function is linear, so we know that it is a minimal surface
  auto exact = [](Eigen::Vector2d x) -> double { return 3 * x[0] - 2 * x[1]; };

  // Compute the minimal surface
  Eigen::VectorXd mu_vec = MinimalGraphSurface::graphMinimalSurface(
      fe_space, exact, 1e-7, 1e-7, 100);

  // We generate global mesh functions to compute the error
  lf::mesh::utils::MeshFunctionGlobal exact_mf{exact};
  const lf::fe::MeshFunctionFE mu_fe(fe_space, mu_vec);

  // The computed minimal surface should agree with the actual minimal surface
  double L2err = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(exact_mf - mu_fe), 2));
  double tol = 1e-8;
  EXPECT_NEAR(L2err, 0., tol);
}

TEST(MinimalGraphSurface, graphMinimalSurface2) {
  // We generate a mesh of a simple 3x3 domain
  auto mesh_0_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // The generated mesh is too course . Therefore, we refine the mesh.
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_0_p, 6);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Generate non-linear boundary data
  auto boundary_data = [](Eigen::Vector2d x) -> double {
    return x[0] * x[0] * x[1] + std::cos(x[0] * x[1]) * cos(x[1] * x[1]);
  };

  // Compute the minimal surface
  Eigen::VectorXd mu_vec = MinimalGraphSurface::graphMinimalSurface(
      fe_space, boundary_data, 1e-12, 1e-12, 100);

  // Construct the DoFHandler to loop over the DoFs on the vertices
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // We will make perturbations to the minimal surface only in the interior of
  // the domain
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 2)};

  // Compute the minimal surface
  double min_surface = MinimalGraphSurface::computeGraphArea(fe_space, mu_vec);

  // We compute several perturbations to the interior of the domain and compute
  // the area of the generated surface. This area should be bigger than the area
  // of the minimal surface.
  for (int i = 0; i < 100; ++i) {
    // Generate random perturbations
    Eigen::VectorXd pert = .1 * Eigen::VectorXd::Random(mu_vec.size());

    // Make sure there are no perturbations on the boundary
    for (int j = 0; j < pert.size(); ++j)
      if (bd_flags(dofh.Entity(j))) pert(j) = 0.;

    // Add the perturbations to the minimal surface
    Eigen::VectorXd mu_vec_perturbed = mu_vec + pert;

    // Compute the area of the perturbed surface
    double surface =
        MinimalGraphSurface::computeGraphArea(fe_space, mu_vec_perturbed);

    // Check that the area of the perturbed surface is bigger than the area of
    // the minimal surface.
    EXPECT_GT(surface, min_surface);
  }
}

}  // namespace MinimalGraphSurface::test
