/**
 * @file nonlinschroedingerequation_test.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "../nonlinschroedingerequation.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <complex>

#include "../propagator.h"

namespace NonLinSchroedingerEquation::test {

constexpr std::complex<double> _i(0, 1);

Eigen::VectorXcd create_mu() {
  int N_dofs = 13;

  Eigen::VectorXd mu_real(N_dofs);
  mu_real << -4.1496914847821, 0.192858374584309, 3.56469932115993,
      -3.18434032487816, 1.27079360432336, 0.398221131357874,
      -0.0149104581241818, -0.486016878008117, 3.27056201654218,
      -4.03999102814268, 4.13035390522648, -0.382981058997326,
      -3.13516563507752;

  Eigen::VectorXd mu_imag(N_dofs);
  mu_imag << 0.638829423213992, 0.0352571144640403, -0.863467224883632,
      2.10972704390233, -3.56923020597589, 4.45620403652878,
      -0.0119650049478216, 0.824035557702947, -1.933150237646, -3.5424244939791,
      -0.561002581688212, 0.16715551617804, 0.1712965388913;

  Eigen::VectorXcd mu = mu_real + _i * mu_imag;

  return mu;
}

Eigen::SparseMatrix<double> create_D() {
  int N_dofs = 13;

  Eigen::SparseMatrix<double> D(N_dofs, N_dofs);
  D.diagonal() << 0.5, 0.666666666666667, 0.416666666666667, 1.08333333333333,
      0.75, 0.416666666666667, 0.666666666666667, 1.25, 1.04166666666667,
      0.333333333333333, 0.5, 0.833333333333333, 0.541666666666667;

  return D;
}

Eigen::SparseMatrix<double> create_A() {
  // Triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Stiffness matrix
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix stiffness_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffness_emp, A_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();

  return A;
}

TEST(NonLinSchroedingerEquation, MassElementMatrixProvider) {
  // Triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // My real mass matrix
  lf::assemble::COOMatrix<double> D_COO(N_dofs, N_dofs);
  NonLinSchroedingerEquation::MassElementMatrixProvider mass_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_emp, D_COO);
  Eigen::SparseMatrix<double> D = D_COO.makeSparse();

  // Reference real mass matrix
  Eigen::SparseMatrix<double> D_ref = create_D();

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (D - D_ref).squaredNorm(), tol);
}

TEST(NonLinSchroedingerEquation, Norm) {
  Eigen::VectorXcd mu = create_mu();
  Eigen::SparseMatrix<double> D = create_D();

  // My norm
  double norm = Norm(mu, D);

  // Reference norm
  double norm_ref = 9.4516312525052;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, norm - norm_ref, tol);
}

TEST(NonLinSchroedingerEquation, KineticEnergy) {
  Eigen::VectorXcd mu = create_mu();
  Eigen::SparseMatrix<double> A = create_A();

  // My kinetic energy
  double kineticEnergy = KineticEnergy(mu, A);

  // Reference kinetic energy
  double kineticEnergy_ref = 185.85719401566;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, kineticEnergy - kineticEnergy_ref, tol);
}

TEST(NonLinSchroedingerEquation, InteractionEnergy) {
  Eigen::VectorXcd mu = create_mu();
  Eigen::SparseMatrix<double> D = create_D();

  // My interaction energy
  double interactionEnergy = InteractionEnergy(mu, D);

  // Reference interaction energy
  double interactionEnergy_ref = 370.597523896931;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, interactionEnergy - interactionEnergy_ref, tol);
}

TEST(NonLinSchroedingerEquation, KineticPropagator) {
  Eigen::VectorXcd mu0 = create_mu();
  Eigen::SparseMatrix<double> A = create_A();
  Eigen::SparseMatrix<std::complex<double>> M = _i * create_D();
  double tau = 0.5;

  // My propagated mu0
  KineticPropagator kineticPropagator(A, M, tau);
  Eigen::VectorXcd mu1 = kineticPropagator(mu0);

  // Reference propagated mu0
  int N_dofs = mu0.size();
  Eigen::VectorXd mu1_ref_real(N_dofs);
  mu1_ref_real << -3.56947555935825, -0.899481130498679, 1.5037668895824,
      1.52054950607771, -3.53782254622478, 4.69451086154526, 0.0322136990251863,
      0.563187056018217, -1.37681884434939, -2.67718731995605, 3.697139063251,
      1.16409150537193, -3.15837046330101;
  Eigen::VectorXd mu1_ref_imag(N_dofs);
  mu1_ref_imag << 3.38175802732124, -0.750151600629607, -2.53402065979762,
      2.49681367346699, 0.96057685919604, -2.27253221774087, 0.0475834865954048,
      -0.883190016287789, -2.79323056125239, 5.15288941038818,
      -3.28248771823266, -0.0425049970290095, 1.02499143544504;
  Eigen::VectorXcd mu1_ref = mu1_ref_real + _i * mu1_ref_imag;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu1 - mu1_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(NonLinSchroedingerEquation, InteractionPropagator) {
  Eigen::VectorXcd mu0 = create_mu();
  double tau = 1.0;

  // My propagated mu0
  InteractionPropagator interactionPropagator(tau);
  Eigen::VectorXcd mu1 = interactionPropagator(mu0);

  // Reference propagated mu0
  int N_dofs = mu0.size();
  Eigen::VectorXd mu1_ref_real(N_dofs);
  mu1_ref_real << -2.02038010911912, 0.194070782992831, 1.58498317601217,
      3.29216975114174, -3.75919937267198, 4.25397749072451,
      -0.0149148301353888, 0.356950494994908, -2.80431837652695,
      5.3319156923224, 0.95346317103043, -0.348117389186395, 2.77273009303072;
  Eigen::VectorXd mu1_ref_imag(N_dofs);
  mu1_ref_imag << -3.68050358857491, 0.0278199201325495, -3.30764043246613,
      1.93716024997555, -0.471954155127836, 1.38564425197975,
      -0.0119595546282825, 0.887599769147364, -2.56309269992962,
      0.664058621646174, 4.05776481089886, 0.231148743205556, -1.47328676597906;
  Eigen::VectorXcd mu1_ref = mu1_ref_real + _i * mu1_ref_imag;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu1 - mu1_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(NonLinSchroedingerEquation, SplitStepPropagator) {
  Eigen::VectorXcd mu0 = create_mu();
  Eigen::SparseMatrix<double> A = create_A();
  Eigen::SparseMatrix<std::complex<double>> M = _i * create_D();
  double tau = 1.0;

  // My propagated mu0
  SplitStepPropagator splitStepPropagator(A, M, tau);
  Eigen::VectorXcd mu1 = splitStepPropagator(mu0);

  // Reference propagated mu0
  int N_dofs = mu0.size();
  Eigen::VectorXd mu1_ref_real(N_dofs);
  mu1_ref_real << -3.53454531029302, -0.0845223412116841, -2.28539603933269,
      -2.92264362555127, 2.20309190220553, -3.0326683108532, 1.24630722905565,
      -0.118039919478382, 3.63183995368114, -4.71176416835282, 3.84566240825742,
      0.386053731999964, -0.441101360403472;
  Eigen::VectorXd mu1_ref_imag(N_dofs);
  mu1_ref_imag << 2.05069831570019, 0.21046236616285, 1.4768731863614,
      -0.984097910136541, -0.977441622733153, 4.0262428321952,
      0.247091542946164, 0.585793525242502, -1.8658312924731, -6.54809413043552,
      -2.25972859122154, 0.423834193838372, -1.7767992773605;
  Eigen::VectorXcd mu1_ref = mu1_ref_real + _i * mu1_ref_imag;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu1 - mu1_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NonLinSchroedingerEquation::test
