/**
 * @file fisherkpp_test.cc
 * @brief NPDE homework NewProblem code
 * @author Louis Hurschler
 * @date 26.03.2024
 * @copyright Developed at ETH Zurich
 */
// In the interest of not changing the problem text over at
// https://gitlab.math.ethz.ch/ralfh/npdeflipped, I will tell clang-tidy to
// ignore the bugprone-suspicious include warning. (Manuel Saladin, 2024-05-28)
#include "../fisherkpp.cc"  // NOLINT(bugprone-suspicious-include)

#include <Eigen/src/Core/util/Constants.h>
#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>

namespace FisherKPP::test {

TEST(FisherKPP, assembleGalerkinMatrices) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/simple.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Dofhandler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  auto c = [](const Eigen::Vector2d& /*x*/) -> double { return 1.; };
  auto galerkinpair = FisherKPP::assembleGalerkinMatrices(dofh, c);

  Eigen::MatrixXd A(galerkinpair.first);

  ASSERT_TRUE(A.rows() == 5);
  ASSERT_TRUE(A.cols() == 5);

  Eigen::MatrixXd A_ref(5, 5);
  A_ref << 1., 0., 0., 0., -1., 0., 1., 0., 0., -1., 0., 0., 1., 0., -1., 0.,
      0., 0., 1., -1., -1., -1., -1., -1., 4.;

  Eigen::MatrixXd M(galerkinpair.second);

  ASSERT_TRUE(M.rows() == 5);
  ASSERT_TRUE(M.cols() == 5);

  Eigen::MatrixXd M_ref(5, 5);
  const double val1 = 0.083333333;
  const double val2 = 0.020833333;
  const double val3 = 0.041666666;
  const double val4 = 0.166666666;

  M_ref << val1, val2, 0., val2, val3, val2, val1, val2, 0., val3, 0., val2,
      val1, val2, val3, val2, 0., val2, val1, val3, val3, val3, val3, val3,
      val4;

  double tol = 1.0e-4;
  EXPECT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);
  EXPECT_NEAR(0.0, (M - M_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FisherKPP, DiffusionEvolutionOperator) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  double T = 1.;
  unsigned int m = 100;
  double lambda = 0.;
  auto c = [](const Eigen::Vector2d& /*x*/) -> double { return 1.; };
  StrangSplit strang_split(fe_space, T, m, lambda, c);

  double tau = T / m;
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(dofh.NumDofs());
  mu(5) = 100.;

  Eigen::VectorXd evol_op = strang_split.diffusionEvolutionOperator(tau, mu);
  ASSERT_TRUE(evol_op.size() == 13);

  Eigen::VectorXd evol_op_ref(13);
  evol_op_ref << 0.8860176500536046662, -3.9193973404503719138,
      18.224713672358088701, -0.72439296635627326015, 11.870414469155882387,
      39.907495444159948761, -0.0028779875937879687009, 7.3710633484319352249,
      -0.44284315016375641605, 15.299487443318426472, 0.44341096366437005027,
      -1.5989862656933366836, -2.2945437332444056366;

  double tol = 1.0e-4;

  EXPECT_NEAR(0.0, (evol_op - evol_op_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(FisherKPP, Evolution) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  double T = 1.;
  unsigned int m = 100;
  double lambda = 0.;
  auto c = [](const Eigen::Vector2d& /*x*/) -> double { return 1.; };
  StrangSplit strang_split(fe_space, T, m, lambda, c);

  Eigen::VectorXd mu = Eigen::VectorXd::Zero(13);
  mu(5) = 100.;
  Eigen::VectorXd cap = 0.2 * Eigen::VectorXd::Ones(13);
  Eigen::VectorXd evolution_res = strang_split.Evolution(cap, mu);

  Eigen::VectorXd evolution_ref(13);
  evolution_ref << 4.6295029473536501996, 4.629690843145014334,
      4.6298952748606287955, 4.6295700848814229644, 4.6297491853779018101,
      4.629838957610099115, 4.6294508503057025806, 4.62971856057560327,
      4.6295073041784773338, 4.6297816718048032669, 4.6293640663389536982,
      4.6295629146841257295, 4.6297636755536704101;

  double tol = 1.0e-4;
  EXPECT_NEAR(0.0, (evolution_ref - evolution_res).lpNorm<Eigen::Infinity>(),
              tol);
}
}  // namespace FisherKPP::test
