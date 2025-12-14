/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP Test
 * @author Amélie Loher
 * @date 01/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "../coupledsecondorderbvp.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace CoupledSecondOrderBVP::test {

constexpr char mesh_file[] = "meshes/simple.msh";

TEST(CoupledSecondOrderBVP, dropMatrixRowsAndColumns) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space =
      std::make_shared<CoupledSecondOrderBVP::FeSpaceLagrangeO2<double>>(
          mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Obtain an array of boolean flags for the nodes of the mesh, 'true'
  // indicates that the node lies on the boundary
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Similarly for edges
  auto edges_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsAndColumns
  auto bd_selector = [&nodes_bd_flags, &edges_bd_flags,
                      &dofh](unsigned int idx) -> bool {
    if (dofh.Entity(idx).RefEl() == lf::base::RefElType::kPoint) {
      return nodes_bd_flags(dofh.Entity(idx));
    }
    return edges_bd_flags(dofh.Entity(idx));
  };
  // Coefficients
  auto const_one = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  auto const_zero = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  lf::assemble::COOMatrix<double> A0(N_dofs, N_dofs);
  // Computing A0
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_zero)>
      A0_builder(fe_space, const_one, const_zero);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A0_builder, A0);
  // Test dropMatrixRowsAndColumns
  dropMatrixRowsAndColumns(bd_selector, A0);

  const std::vector<Eigen::Triplet<double>> A0_triplets_vec = A0.triplets();

  // Vector storing the nnz entries of A0
  std::vector<double> triplet_values(44);

  for (auto& triplet : A0_triplets_vec) {
    triplet_values.push_back(triplet.value());
  }
  // Compute norm
  double norm = 0.0;
  for (int i = 0; i < triplet_values.size(); ++i) {
    norm += triplet_values[i] * triplet_values[i];
  }
  norm = std::sqrt(norm);

  // Reference values for nnz entries of A0
  Eigen::MatrixXd refval(1, 44);
  refval << 1, -0.666667, -0.666667, -0.666667, 2.66667, -6.10623e-16,
      -0.666667, -6.10623e-16, 2.66667, 1, -0.666667, -0.666667, -0.666667,
      2.66667, -6.10623e-16, -0.666667, -6.10623e-16, 2.66667, 1, -0.666667,
      -0.666667, -0.666667, 2.66667, -6.10623e-16, -0.666667, -6.10623e-16,
      2.66667, 1, -0.666667, -0.666667, -0.666667, 2.66667, -6.10623e-16,
      -0.666667, -6.10623e-16, 2.66667, 1, 1, 1, 1, 1, 1, 1, 1;

  double tol = 1.0e-4;
  ASSERT_NEAR(refval.norm(), norm, tol);
}

TEST(CoupledSecondOrderBVP, dropMatrixRows) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space =
      std::make_shared<CoupledSecondOrderBVP::FeSpaceLagrangeO2<double>>(
          mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Obtain an array of boolean flags for the nodes of the mesh, 'true'
  // indicates that the node lies on the boundary
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Similarly for edges
  auto edges_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsAndColumns
  auto bd_selector = [&nodes_bd_flags, &edges_bd_flags,
                      &dofh](unsigned int idx) -> bool {
    if (dofh.Entity(idx).RefEl() == lf::base::RefElType::kPoint) {
      return nodes_bd_flags(dofh.Entity(idx));
    }
    return edges_bd_flags(dofh.Entity(idx));
  };

  // Coefficients
  auto const_one = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  auto const_zero = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);  // off-diag blocks
  // Computing M: standard mass matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_zero), decltype(const_one)>
      M_builder(fe_space, const_zero, const_one);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_builder, M);

  // Enforce Dirichlet boundary conditions
  dropMatrixRows(bd_selector, M);
  // Vector of Triplets for M
  const std::vector<Eigen::Triplet<double>> M_triplets_vec = M.triplets();

  // Vector storing the nnz entries of A0
  std::vector<double> triplet_values(72);
  for (auto& triplet : M_triplets_vec) {
    triplet_values.push_back(triplet.value());
  }
  // Compute norm
  double norm = 0;
  for (int i = 0; i < triplet_values.size(); ++i) {
    norm += triplet_values[i] * triplet_values[i];
  }
  norm = std::sqrt(norm);

  // Reference values for nnz entries of A0
  Eigen::MatrixXd refval(1, 72);
  refval << -0.00138889, -0.00138889, 0.00833333, -0.00555556, 1.30104e-17,
      -6.07153e-18, -0.00555556, -6.50521e-18, 1.30104e-17, 0.0222222,
      0.0444444, 0.0222222, 3.03577e-18, -0.00555556, -6.07153e-18, 0.0222222,
      0.0222222, 0.0444444, -0.00138889, -0.00138889, 0.00833333, -0.00555556,
      1.30104e-17, -6.07153e-18, -0.00555556, -6.50521e-18, 1.30104e-17,
      0.0222222, 0.0444444, 0.0222222, 3.03577e-18, -0.00555556, -6.07153e-18,
      0.0222222, 0.0222222, 0.0444444, -0.00138889, -0.00138889, 0.00833333,
      -0.00555556, 1.30104e-17, -6.07153e-18, -0.00555556, -6.50521e-18,
      1.30104e-17, 0.0222222, 0.0444444, 0.0222222, 3.03577e-18, -0.00555556,
      -6.07153e-18, 0.0222222, 0.0222222, 0.0444444, -0.00138889, -0.00138889,
      0.00833333, -0.00555556, 1.30104e-17, -6.07153e-18, -0.00555556,
      -6.50521e-18, 1.30104e-17, 0.0222222, 0.0444444, 0.0222222, 3.03577e-18,
      -0.00555556, -6.07153e-18, 0.0222222, 0.0222222, 0.0444444;

  double tol = 1.0e-4;
  ASSERT_NEAR(refval.norm(), norm, tol);
}

TEST(CoupledSecondOrderBVP, solveCoupledBVP) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space =
      std::make_shared<CoupledSecondOrderBVP::FeSpaceLagrangeO2<double>>(
          mesh_p);
  // reaction coefficient
  double gamma = 1.0;
  // Right-hand side source function f
  auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::cos(x.norm()); });

  // Solution Vector
  Eigen::VectorXd sol = solveCoupledBVP(fe_space, gamma, f);

  // Reference Solution
  Eigen::VectorXd ref_sol(26);
  ref_sol << 0, 0, 0, 0, 0.0472408, 0.0712641, 0.0417256, 0.0417256, 0.0712641,
      0.0402412, 0.0344878, 0.0290578, 0.0344878, -0.0227536, -0.022412,
      -0.0220806, -0.022412, -0.0237479, -0.0230559, -0.0227161, -0.0227161,
      -0.0230559, -0.0233836, -0.0231244, -0.0228704, -0.0231244;

  double tol = 1.0e-6;
  ASSERT_NEAR(ref_sol.norm(), sol.norm(), tol);
}

}  // namespace CoupledSecondOrderBVP::test
