/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <memory>

#include "../../meshes/mesh.h"
#include "../mylinearfeelementmatrix.h"
#include "../mylinearloadvector.h"
#include "../solve.h"

namespace ElementMatrixComputation::test {

double numericalPrecision = 1e-8;

double f(Eigen::Vector2d x) { return 1 + x(0) * x(0) + x(1) * x(1); };

// useful functors
Eigen::Matrix<double, 2, 2> zeroMatrixFunctor(Eigen::Vector2d x) {
  return Eigen::MatrixXd::Zero(2, 2);
};

double zeroScalar(Eigen::Vector2d x) { return 0.0; };

Eigen::Matrix<double, 2, 2> identityMatrixFunctor(Eigen::Vector2d x) {
  return Eigen::MatrixXd::Identity(2, 2);
};

double identityScalarFunctor(Eigen::Vector2d x) { return 1.0; };

//////////////////////
// Test solve
//////////////////////
TEST(Solve, test) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(8, 1.0 / 3.0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs(dofh.NumDofs());

  lf::mesh::utils::MeshFunctionGlobal mf_alpha{identityMatrixFunctor};
  lf::mesh::utils::MeshFunctionGlobal mf_gamma{identityScalarFunctor};

  // Galerkin matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_alpha), decltype(mf_gamma)>
      elmat_builder(fe_space, mf_alpha, mf_gamma);

  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Right hand side
  lf::mesh::utils::MeshFunctionGlobal mf_f{ElementMatrixComputation::f};

  lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_f)> elvec_builder(
      mf_f);
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // solve
  Eigen::VectorXd sol_vec = Eigen::VectorXd::Zero(N_dofs);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  sol_vec = solver.solve(phi);

  double error =
      (sol_vec - ElementMatrixComputation::solve(elmat_builder, elvec_builder))
          .norm();

  EXPECT_LT(error, numericalPrecision);
}

//////////////////////
// Test SolvePoissonBVP
//////////////////////
TEST(SolvePoissonBVP, test) {
  lf::mesh::utils::MeshFunctionGlobal mf_f{ElementMatrixComputation::f};

  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder;
  lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_f)> elvec_builder(
      mf_f);

  Eigen::VectorXd solution =
      ElementMatrixComputation::solve(elmat_builder, elvec_builder);

  Eigen::VectorXd student_solution =
      ElementMatrixComputation::solvePoissonBVP();

  if (student_solution.norm() < numericalPrecision) {
    EXPECT_TRUE(false);
  } else {
    double error = (solution - student_solution).norm();

    EXPECT_LT(error, numericalPrecision);
  }
}

//////////////////////
// Test MyLinearLoadVector
//////////////////////
TEST(MyLinearLoadVector, testTriangles) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);

  lf::mesh::utils::MeshFunctionGlobal mf_f{ElementMatrixComputation::f};

  ElementMatrixComputation::MyLinearLoadVector elvec_builder(
      ElementMatrixComputation::f);
  lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_f)>
      elvec_builder_exact(mf_f);

  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    if (cell->RefEl() == lf::base::RefEl::kTria()) {
      auto elem_vec = elvec_builder.Eval(*cell);
      auto elem_vec_exact = elvec_builder_exact.Eval(*cell);
      double error = (elem_vec - elem_vec_exact).norm();
      EXPECT_LT(error, numericalPrecision);
    }
  }
}

TEST(MyLinearLoadVector, testQuads) {
  auto mesh_p = Generate2DTestMesh();

  lf::mesh::utils::MeshFunctionGlobal mf_f{ElementMatrixComputation::f};

  ElementMatrixComputation::MyLinearLoadVector elvec_builder(
      ElementMatrixComputation::f);
  lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_f)>
      elvec_builder_exact(mf_f);

  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    if (cell->RefEl() == lf::base::RefEl::kQuad()) {
      auto elem_vec = elvec_builder.Eval(*cell);
      auto elem_vec_exact = elvec_builder_exact.Eval(*cell);
      double error = (elem_vec - elem_vec_exact).norm();
      EXPECT_LT(error, numericalPrecision);
    }
  }
}

//////////////////////
// Test MyLinearFEElementMatrix
//////////////////////
TEST(MyLinearFEElementMatrix, testTriangles) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};

  lf::mesh::utils::MeshFunctionGlobal mf_alpha{identityMatrixFunctor};
  lf::mesh::utils::MeshFunctionGlobal mf_gamma{identityScalarFunctor};

  ElementMatrixComputation::MyLinearFEElementMatrix elmat_builder;
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_alpha), decltype(mf_gamma)>
      elmat_builder_exact(fe_space, mf_alpha, mf_gamma);

  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    if (cell->RefEl() == lf::base::RefEl::kTria()) {
      auto elem_mat = elmat_builder.Eval(*cell);
      auto elem_mat_exact = elmat_builder_exact.Eval(*cell);
      double error = (elem_mat.block(0, 0, 3, 3) - elem_mat_exact).norm();
      EXPECT_LT(error, numericalPrecision);
    }
  }
}

TEST(MyLinearFEElementMatrix, testQuads) {
  // std::shared_ptr<lf::mesh::Mesh> mesh_p =
  //     lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  auto mesh_p = Generate2DTestMesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};

  lf::mesh::utils::MeshFunctionGlobal mf_alpha{identityMatrixFunctor};
  lf::mesh::utils::MeshFunctionGlobal mf_gamma{identityScalarFunctor};

  ElementMatrixComputation::MyLinearFEElementMatrix elmat_builder;
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_alpha), decltype(mf_gamma)>
      elmat_builder_exact(fe_space, mf_alpha, mf_gamma);

  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    if (cell->RefEl() == lf::base::RefEl::kQuad()) {
      auto elem_mat = elmat_builder.Eval(*cell);
      auto elem_mat_exact = elmat_builder_exact.Eval(*cell);
      double error = (elem_mat - elem_mat_exact).norm();
      EXPECT_LT(error, numericalPrecision);
    }
  }
}

//////////////////////
// Test solveNeumannEq
//////////////////////
TEST(solveNeumannEq, test) {
  ElementMatrixComputation::MyLinearFEElementMatrix elmat_builder;
  ElementMatrixComputation::MyLinearLoadVector elvec_builder(
      ElementMatrixComputation::f);

  Eigen::VectorXd solution =
      ElementMatrixComputation::solve(elmat_builder, elvec_builder);
  Eigen::VectorXd student_solution = ElementMatrixComputation::solveNeumannEq();

  if (student_solution.norm() < numericalPrecision) {
    EXPECT_TRUE(false);
  } else {
    double error = (solution - student_solution).norm();

    EXPECT_LT(error, numericalPrecision);
  }
}

}  // namespace ElementMatrixComputation::test
