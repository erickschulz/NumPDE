/** @file
 * @brief NPDE SDIRKMethodOfLines Test
 * @author Amélie Loher
 * @date 07/04/2020
 * @copyright Developed at ETH Zurich
 */

#include "../sdirkmethodoflines.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

// In the interest of not changing the problem text over at
// https://gitlab.math.ethz.ch/ralfh/npdeflipped, I will tell clang-tidy to
// ignore the bugprone-suspicious-include warning. (Manuel Saladin, 2024-05-28)
#include "../sdirkmethodoflines_ode.cc"  // NOLINT(bugprone-suspicious-include)

namespace SDIRKMethodOfLines::test {

TEST(SDIRKMethodOfLines, assembleGalerkinMatrices) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  double c = 1.0;

  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> pair =
      assembleGalerkinMatrices(dofh, c);

  Eigen::MatrixXd A(pair.first);
  Eigen::MatrixXd M(pair.second);

  Eigen::MatrixXd A_ref(N_dofs, N_dofs);

  A_ref << 1.66667, 0.166667, 0, 0.166667, -1, 0.166667, 1.66667, 0.166667, 0,
      -1, 0, 0.166667, 1.66667, 0.166667, -1, 0.166667, 0, 0.166667, 1.66667,
      -1, -1, -1, -1, -1, 4;

  Eigen::MatrixXd M_ref(N_dofs, N_dofs);

  M_ref << 0.0833333, 0.0208333, 0, 0.0208333, 0.0416667, 0.0208333, 0.0833333,
      0.0208333, 0, 0.0416667, 0, 0.0208333, 0.0833333, 0.0208333, 0.0416667,
      0.0208333, 0, 0.0208333, 0.0833333, 0.0416667, 0.0416667, 0.0416667,
      0.0416667, 0.0416667, 0.16667;

  double tol = 1.0e-4;

  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (M - M_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(SDIRKMethodOfLines, solveTemperatureEvolution) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  unsigned int m = 12;

  double c = 1.0;

  Eigen::VectorXd init(N_dofs);
  for (int idx = 0; idx < N_dofs; idx++) {
    init(idx) = 5.0;
  }

  std::pair<Eigen::VectorXd, Eigen::VectorXd> pair =
      solveTemperatureEvolution(dofh, m, c, init);

  Eigen::VectorXd sol = pair.first;
  Eigen::VectorXd eng = pair.second;

  Eigen::VectorXd sol_ref(N_dofs);
  Eigen::VectorXd eng_ref(m + 1);

  sol_ref << 0.122308, 0.122308, 0.122308, 0.122308, 0.165211;
  eng_ref << 5, 3.68005, 2.73267, 2.02427, 1.50051, 1.11206, 0.824213, 0.610865,
      0.452744, 0.335552, 0.248695, 0.18432, 0.136609;

  double tol = 1.0e-4;

  ASSERT_NEAR(0.0, (sol - sol_ref).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (eng - eng_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(SDIRKMethodOfLines, thermalEnergy) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  Eigen::VectorXd init(N_dofs);
  for (int idx = 0; idx < N_dofs; idx++) {
    init(idx) = 5.0;
  }

  double thrmeng = thermalEnergy(dofh, init);
  double thrmeng_ref = 5.0;

  double tol = 1.0e-6;

  ASSERT_NEAR(0.0, thrmeng - thrmeng_ref, tol);
}

TEST(SDIRKMethodOfLines, sdirk2SteppingLinScalODE) {
  unsigned int m = 12;

  std::vector<double> sol = sdirk2SteppingLinScalODE(m);

  std::vector<double> sol_ref;

  sol_ref.push_back(1);
  sol_ref.push_back(0.84632);
  sol_ref.push_back(0.716258);
  sol_ref.push_back(0.606184);
  sol_ref.push_back(0.513026);
  sol_ref.push_back(0.434184);
  sol_ref.push_back(0.367459);
  sol_ref.push_back(0.310988);
  sol_ref.push_back(0.263196);
  sol_ref.push_back(0.222748);
  sol_ref.push_back(0.188516);
  sol_ref.push_back(0.159545);
  sol_ref.push_back(0.135026);

  double norm = 0.0;
  double norm_ref = 0.0;
  for (int i = 0; i < sol.size(); i++) {
    norm += sol[i] * sol[i];
    norm_ref += sol_ref[i] * sol_ref[i];
  }

  norm = std::sqrt(norm);
  norm_ref = std::sqrt(norm_ref);

  double tol = 1.0e-4;
  ASSERT_NEAR(0.0, norm - norm_ref, tol);
}

} /* namespace SDIRKMethodOfLines::test */
