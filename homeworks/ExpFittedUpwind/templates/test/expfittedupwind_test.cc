/**
 * @file expfittedupwind_test.cc
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher, Philippe Peter
 * @date 07.01.2021
 * @copyright Developed at ETH Zurich
 */

#include "../expfittedupwind.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <memory>
#include <vector>

namespace ExpFittedUpwind::test {

// The Bernoulli function should be monotonically decreasing
TEST(Bernoulli, Monotinicity) {
  std::vector<double> values = {3.0,     1.5,    0.5,     1.0E-5,
                                1.0E-7,  1.0E-9, -1.0E-9, -1.0E-7,
                                -1.0E-5, -1.0,   -5.0,    -6.0};

  for (int i = 0; i < values.size() - 1; ++i) {
    EXPECT_LE(Bernoulli(values[i]), Bernoulli(values[i + 1]));
  }
}

// Check that B(tau) does not explode for small tau
TEST(Bernoulli, Cancellation) {
  for (double i = 8.0; i >= std::numeric_limits<double>::min(); i /= 2.0) {
    EXPECT_LE(std::abs(Bernoulli(i)), 2.0);
  }
}

// Check that B(log(2)) = log(2)
TEST(Bernoulli, Evaluation) {
  EXPECT_NEAR(Bernoulli(std::log(2)), std::log(2), 1.0E-8);
}

// For Psi = c, all entries of beta are given by beta(e) = std::exp(c)
TEST(CompBeta, ConstantPSI) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d /*x*/) { return 3.1; });
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);
  double ref = std::exp(3.1);

  auto beta = CompBeta(*fe_space, mu);

  for (auto entity : mesh_p->Entities(1)) {
    EXPECT_DOUBLE_EQ((*beta)(*entity), ref);
  }
}

// Verify the computation of beta for a linear psi
// based on precomputed reference values on a test mesh.
TEST(CompBeta, linearPSI) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return 1.0 + x(0) + 2 * x(1); });
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);

  auto beta = CompBeta(*fe_space, mu);
  auto edges = mesh_p->Entities(1);

  // edges of triangle 0
  EXPECT_DOUBLE_EQ((*beta)(*(edges[0])), std::exp(2.5) * Bernoulli(1.5));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[1])), std::exp(4.0) * Bernoulli(3.0));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[2])), std::exp(4.0) * Bernoulli(3.0));

  // edges of triangle 6:
  EXPECT_DOUBLE_EQ((*beta)(*(edges[12])), std::exp(6.0) * Bernoulli(1.0));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[13])), std::exp(6.0) * Bernoulli(1.0));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[14])), std::exp(6.0) * Bernoulli(0.0));

  // edges of triangle 11
  EXPECT_DOUBLE_EQ((*beta)(*(edges[18])), std::exp(6.0) * Bernoulli(0.0));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[20])), std::exp(6.0) * Bernoulli(-2.5));
  EXPECT_DOUBLE_EQ((*beta)(*(edges[23])), std::exp(8.5) * Bernoulli(2.5));
}

// for Psi = const = c, the exponentially fitted element matrix reduces to
// A^{exp}_K = -A_K where A_K is the element matrix for (u,v) -> \int_K grad u *
// grad v dx
TEST(ExpFittedEMP, Psi_const) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return 3.0; });
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);

  // Implemented exponentially fitted upwind provider
  ExpFittedEMP upwind_provider(fe_space, mu);

  // Lehrfem++ element matrix provider for the the minus laplacian
  lf::uscalfe::LinearFELaplaceElementMatrix standard_provider;

  // Expect that for Psi=const the ExpFitted Element matrix is A_k
  for (auto* cell : mesh_p->Entities(0)) {
    Eigen::Matrix3d E_k = upwind_provider.Eval(*cell);
    Eigen::Matrix3d A_k = standard_provider.Eval(*cell).block<3, 3>(0, 0);

    EXPECT_NEAR((E_k - A_k).norm(), 0.0, 1.0E-10);
  }
}

// The last two tests are based on the fact, that
// the expontially Fitted upwind scheme provides a system matrix  A that
// corresponds to the bilinear form (u,v) -> b(u,v) \int_{\Omega} -j(u,Psi) *
// grad v dx (j(u,Psi) = const ) In particular this value can be approximated by
// b_v^T * A * b_u, where b_v and b_u are the nodal projections of v and u into
// the underlying FE space

// u = exp(Psi)
// j(u,Psi) = 0
// and b(u,v) = 0 for any v
TEST(ExpFittedEMP, Bilinear_form_1) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  Eigen::Vector2d q = Eigen::Vector2d::Ones(2);
  auto Psi = [&q](Eigen::Vector2d x) { return q.dot(x); };

  auto u = [&Psi](Eigen::Vector2d x) { return std::exp(Psi(x)); };

  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(Psi);
  auto mf_u = lf::mesh::utils::MeshFunctionGlobal(u);

  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);
  Eigen::VectorXd u_vec = lf::fe::NodalProjection(*fe_space, mf_u);

  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  ExpFittedEMP elmat_builder(fe_space, mu);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  EXPECT_NEAR((A_crs * u_vec).norm(), 0.0, 1.0E-10);
}

// Psi = q' * x
// u = c (--> j(u) = c*q)
// v = r' *x
// b(u,v) = -|\Omega|*c*q'*r
TEST(ExpFittedEMP, Bilinear_form_2) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  double area = 9.0;
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  Eigen::Vector2d q = Eigen::Vector2d::Ones(2);
  auto Psi = [&q](Eigen::Vector2d x) { return q.dot(x); };
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(Psi);

  double u_value = 2.0;
  auto mf_u = lf::mesh::utils::MeshFunctionConstant(u_value);

  Eigen::Vector2d r = Eigen::Vector2d::Ones(2);
  auto v = [&r](Eigen::Vector2d x) { return r.dot(x); };
  auto mf_v = lf::mesh::utils::MeshFunctionGlobal(v);

  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);
  Eigen::VectorXd u_vec = lf::fe::NodalProjection(*fe_space, mf_u);
  Eigen::VectorXd v_vec = lf::fe::NodalProjection(*fe_space, mf_v);

  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  ExpFittedEMP elmat_builder(fe_space, mu);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  double value = v_vec.transpose() * A_crs * u_vec;
  double reference = -area * u_value * r.dot(q);
  EXPECT_NEAR(value, reference, 1.0E-10);
}
// Psi = q * x
// dirichlet bc with g = exp(Psi(x))
// => u(x) = exp(Psi(x))
// Note that this test would not detect a sign error in the
// eval function of the ExpFittedEMP class because f == 0 and we provide proper
// dirichlet data
TEST(ExpFittedEMP, SolveDriftEquation) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::Vector2d q = Eigen::Vector2d::Ones(2);
  auto Psi = [&q](const Eigen::Vector2d& x) { return q.dot(x); };
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(Psi);
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);

  auto exp_psi = [&Psi](const Eigen::VectorXd& x) { return std::exp(Psi(x)); };
  auto f_zero = [](const Eigen::VectorXd& /*x*/) { return 0.; };

  Eigen::VectorXd res =
      solveDriftDiffusionDirBVP(fe_space, mu, f_zero, exp_psi);

  auto mf_exp_psi = lf::mesh::utils::MeshFunctionGlobal(exp_psi);
  Eigen::VectorXd res_analytically =
      lf::fe::NodalProjection(*fe_space, mf_exp_psi);

  EXPECT_NEAR((res - res_analytically).norm(), 0.0, 1e-10);
}

// Test with Phi == const, f == c
// In this case, we solve the equation
// div(- grad u) = c => u = -(c * (x_1^2 + x_2^2))/4.
// use g = u on boundary
TEST(ExpFittedEMP, TestPhiZero) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/island.msh");

  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto Psi = [](const Eigen::Vector2d& /*x*/) { return 0.; };
  auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal(Psi);
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space, mf_Psi);

  const double c = 3.;
  auto f_const = [c](const Eigen::VectorXd& /*x*/) { return c; };

  auto g = [&c](const Eigen::Vector2d& x) { return -c * x.dot(x) * .25; };

  Eigen::VectorXd res = solveDriftDiffusionDirBVP(fe_space, mu, f_const, g);

  auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);
  Eigen::VectorXd res_analytically = lf::fe::NodalProjection(*fe_space, mf_g);

  // note that we have to choose a high tolerance here because the
  // solution is quadratic and we interpolate it with piecewise linear funcitons
  // on a relatively big mesh
  EXPECT_NEAR((res - res_analytically).lpNorm<Eigen::Infinity>(), 0.0, 1e-1);
}

}  // namespace ExpFittedUpwind::test
