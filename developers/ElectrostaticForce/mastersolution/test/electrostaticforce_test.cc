/**
 * @file electrostaticforce.cc
 * @brief ElectrostaticForce
 * @author Erick Schulz
 * @date 27.11.2019
 * @copyright Developed at ETH Zurich
 */

// HACK:
#undef SOLUTION
#define SOLUTION 1

#include "../electrostaticforce.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace ElectrostaticForce::test {

TEST(ElectrostaticForce, computeExactForce) {
  Eigen::Vector2d exactForce = ElectrostaticForce::computeExactForce();
  double tol = 1.0e-3;
  ASSERT_NEAR(13.0776, exactForce(0), tol);
  ASSERT_NEAR(0.0, exactForce(1), tol);
}

TEST(ElectrostaticForce, computeForceDomainFunctional) {
  std::string mesh_file = "meshes/emforce" + std::to_string(4) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::VectorXd approx_sol = ElectrostaticForce::solvePoissonBVP(fe_space_p);

  Eigen::Vector2d approx_force_domain_functional =
      ElectrostaticForce::computeForceDomainFunctional(fe_space_p, approx_sol);

  double tol = 1.0e-4;
  ASSERT_NEAR(13.0711, approx_force_domain_functional(0), tol);
  ASSERT_NEAR(0.01106, approx_force_domain_functional(1), tol);
}

TEST(ElectrostaticForce, computeForceBoundaryFunctional) {
  std::string mesh_file = "meshes/emforce" + std::to_string(4) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::VectorXd approx_sol = ElectrostaticForce::solvePoissonBVP(fe_space_p);

  Eigen::Vector2d approx_force_boundary_functional =
      ElectrostaticForce::computeForceBoundaryFunctional(fe_space_p,
                                                         approx_sol);

  double tol = 1.0e-3;
  ASSERT_NEAR(12.5894, approx_force_boundary_functional(0), tol);
  ASSERT_NEAR(-0.0602, approx_force_boundary_functional(1), tol);
}

TEST(ElectrostaticForce, solvePoissonBVPBoundaryConditions) {
  double tol = 1.0e-14;
  std::string mesh_file = "meshes/emforce" + std::to_string(1) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};

  Eigen::VectorXd approx_sol = ElectrostaticForce::solvePoissonBVP(fe_space_p);

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    if (bd_flags(*node)) {
      auto dof_idx = dofh.GlobalDofIndices(*node);
      auto endpoints = lf::geometry::Corners(*(node->Geometry()));
      if (endpoints.norm() < 0.27) {
        ASSERT_NEAR(1.0, approx_sol(dof_idx[0]), tol);
      } else {
        ASSERT_NEAR(0.0, approx_sol(dof_idx[0]), tol);
      }
    }
  }
}

}  // namespace ElectrostaticForce::test
