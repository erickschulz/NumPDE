/**
 * @file solavgboundary_test.cc
 * @brief NPDE homework Solavgboundary code
 * @author Bob Schreiner
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../solavgboundary.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace solavgboundary::test {

constexpr char mesh_file[] = "meshes/square.msh";

constexpr auto const_one = [](Eigen::Vector2d x) -> double { return 1.0; };
constexpr auto const_zero = [](Eigen::Vector2d x) -> double { return 0.0; };

TEST(AugmentedMatrixTest, augmentMatrix) {
  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  const lf::assemble::DofHandler* dofh = &fe_space->LocGlobMap();
  const int N_dofs = dofh->NumDofs();

  // Compute c explicitly and not rely on an implemented function
  const Eigen::SparseMatrix<double> B =
      compGalerkinMatrix(*dofh, const_zero, const_zero, const_one).makeSparse();
  const Eigen::VectorXd c = B * Eigen::VectorXd::Constant(N_dofs, 1.);
  const Eigen::SparseMatrix<double> A_test = augmentMatrix(fe_space, c);

  // Compute A with the mastersolution
  lf::assemble::COOMatrix<double> A =
      lf::assemble::COOMatrix<double>(N_dofs + 1, N_dofs + 1);

  A = compGalerkinMatrix(*dofh, const_one, const_zero, const_zero);
  for (int i = 0; i < N_dofs; ++i) {
    A.AddToEntry(N_dofs, i, c[i]);
    A.AddToEntry(i, N_dofs, c[i]);
  }
  const Eigen::SparseMatrix<double> A_sprs = A.makeSparse();

  ASSERT_NEAR(A_sprs.norm(), A_test.norm(), 1e-8);
}

TEST(VectorCTest, computeCVector) {
  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  const Eigen::VectorXd c_test = computeCVector(fe_space);

  const lf::assemble::DofHandler* dofh = &fe_space->LocGlobMap();
  const int N_dofs = dofh->NumDofs();

  const Eigen::SparseMatrix<double> B =
      compGalerkinMatrix(*dofh, const_zero, const_zero, const_one).makeSparse();
  const Eigen::VectorXd c = B * Eigen::VectorXd::Constant(N_dofs, 1.);
  ASSERT_NEAR(c.norm(), c_test.norm(), 1e-8);
}
}  // namespace solavgboundary::test
