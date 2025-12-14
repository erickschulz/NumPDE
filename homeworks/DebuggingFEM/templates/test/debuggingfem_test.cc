#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <functional>

#include "../locallaplaceqfe.h"
#include "../qfeinterpolator.h"
#include "../qfeprovidertester.h"

namespace DebuggingFEM::test {

constexpr double Square(double x) { return x * x; }

struct TestPair {
  typedef Eigen::Matrix<double, 6, 1> Vector6d;

  TestPair(const Vector6d& a) {
    function = [a](Eigen::Vector2d x) {
      return a(0) * Square(x(0)) + a(1) * Square(x(1)) + a(2) * x(0) * x(1) +
             a(3) * x(0) + a(4) * x(1) + a(5);
    };
    // coefficients of quadratic polynomial |grad function|^2
    Vector6d b;
    b << 4.0 * Square(a(0)) + Square(a(2)), 4.0 * Square(a(1)) + Square(a(2)),
        4.0 * a(2) * (a(0) + a(1)), 2.0 * (2.0 * a(0) * a(3) + a(2) * a(4)),
        2.0 * (2.0 * a(1) * a(4) + a(2) * a(3)), Square(a(3)) + Square(a(4));
    // monomials integrated over [0, 1]^2
    Eigen::VectorXd integrated_monomials(6);
    integrated_monomials << 1.0 / 3.0, 1.0 / 3.0, 0.25, 0.5, 0.5, 1.0;
    // reference value for energy (H_1 seminorm)
    energy = integrated_monomials.dot(b);
  }

  std::function<double(Eigen::Vector2d x)> function;
  double energy;
};

TEST(DebuggingFEM, interpolateOntoQuadFE) {
  // reference
  Eigen::VectorXd a(6);
  a << 1.0, 4.0, 2.0, 3.0, 2.0, 1.0;
  TestPair test_pair(a);
  auto p = test_pair.function;
  double energy_ref = test_pair.energy;

  // to test
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1},
                                          {lf::base::RefEl::kSegment(), 1},
                                          {lf::base::RefEl::kTria(), 0},
                                          {lf::base::RefEl::kQuad(), 1}});

  const lf::base::size_type N_dofs(dofh.NumDofs());
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  auto element_matrix_provider = DebuggingFEM::LocalLaplaceQFE2();
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      mat);

  const Eigen::SparseMatrix<double> A(mat.makeSparse());
  Eigen::VectorXd eta = DebuggingFEM::interpolateOntoQuadFE(dofh, p);
  double energy = eta.dot(A * eta);

  EXPECT_NEAR(energy, energy_ref, 1.0e-8);
}

TEST(DebuggingFEM, QFEProviderTester) {
  // reference
  Eigen::VectorXd a(6);
  a << 1.0, 4.0, 2.0, 3.0, 2.0, 1.0;
  TestPair test_pair(a);
  auto p = test_pair.function;
  double energy_ref = test_pair.energy;

  // to test
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1},
                                          {lf::base::RefEl::kSegment(), 1},
                                          {lf::base::RefEl::kTria(), 0},
                                          {lf::base::RefEl::kQuad(), 1}});

  auto element_matrix_provider = DebuggingFEM::LocalLaplaceQFE2();
  QFEProviderTester qfe_provider_tester(dofh, element_matrix_provider);
  double energy = qfe_provider_tester.energyOfInterpolant(p);

  EXPECT_NEAR(energy, energy_ref, 1.0e-8);
}

}  // namespace DebuggingFEM::test
