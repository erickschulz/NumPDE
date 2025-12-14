/**
 * @file LinHypDampWaveSys_test.cc
 * @brief NPDE homework LinHypDampWaveSys code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../linhypdampwavesys.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace LinHypDampWaveSys::test {

TEST(LinHypDampWaveSys, linsysoderhs) {
  unsigned n = 5;                                     // length of state vector
  Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(2, n);   // Discrete state
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2, n);  // Computed RHS
  Eigen::MatrixXd rhs_exact = Eigen::MatrixXd::Zero(2, n);  // exact RHS

  // Initialize discrete state and exact RHS for comparison
  mu << 1, 1, 1, 0, 0, 1, 1, 1, 0, 0;
  rhs_exact << 1., 1., 3.5, -1.5, 1., 0., 0., 2.5, -2.5, 0.;

  // We come up with a flux. Note that this is not a realistic flux
  auto numflux = [](Eigen::VectorXd mu1,
                    Eigen::VectorXd mu2) -> Eigen::VectorXd {
    return -.5 * mu1 + .5 * mu2;
  };

  // We have the following source term
  auto g = [](Eigen::VectorXd u) -> Eigen::VectorXd {
    Eigen::VectorXd out(2);
    out << 1., 0.;
    return out;
  };

  // Compute the RHS
  rhs = linsysoderhs(mu, .2, numflux, g);

  // Make sure the RHS was computed correctly
  ASSERT_LE((rhs - rhs_exact).norm(), 1e-6);
}

TEST(LinHypDampWaveSys, fvEulLinSys) {
  // Some parameters
  double a = -1.;
  double b = 2.;
  double T = 3.;
  unsigned int N = 5;
  unsigned int M = 3;

  // We come up with a flux. Note that this is not a realistic flux
  auto numflux = [](Eigen::VectorXd mu1,
                    Eigen::VectorXd mu2) -> Eigen::VectorXd {
    return -.5 * mu1 + .5 * mu2;
  };

  // We have the following source term
  auto g = [](Eigen::VectorXd u) -> Eigen::VectorXd {
    Eigen::VectorXd out(2);
    out << 1., 0.;
    return out;
  };

  // Define an initial condition
  auto u0 = [](double x) -> Eigen::Vector2d {
    Eigen::Vector2d out;
    out(0) = 0.;
    if (0 <= x && x <= .5)
      out(1) = 1.;
    else if (0.5 <= x && x <= 1)
      out(1) = -1.;
    else
      out(1) = 0.;
    return out;
  };

  // Compute the discrete solution
  auto sol = fvEulLinSys(a, b, T, N, M, u0, numflux, g);

  // The exact solution is given as follows
  Eigen::MatrixXd sol_exact = Eigen::MatrixXd(2, N);
  sol_exact << 3., 3., 3., 3., 3., 4.97685, -19.5139, 30.0741, -19.5139,
      4.97685;

  // Compare exact and computed solutions
  ASSERT_LE((sol - sol_exact).norm(), 1e-3);
}

TEST(LinHypDampWaveSys, ev1DampWave) {
  // Some parameters
  double c = 2.;
  double r = 3.;
  double T = 4.;
  unsigned int N = 5;

  // Define an initial condition
  auto u0 = [](double x) -> Eigen::Vector2d {
    Eigen::Vector2d out;
    out(0) = 0.;
    if (0 <= x && x <= .5)
      out(1) = 1.;
    else if (0.5 <= x && x <= 1)
      out(1) = -1.;
    else
      out(1) = 0.;
    return out;
  };

  // Compute the solution using upwind and Lax-Friedrichs flux
  auto sol_LF = ev1DampWave(c, r, T, N, u0, LF);
  auto sol_UW = ev1DampWave(c, r, T, N, u0, UW);

  // Define the exact solutions
  Eigen::MatrixXd sol_LF_exact(2, N);
  sol_LF_exact << -3.29285, 13.7534, -21.8933, 70.9194, -85.7135, -0.41613,
      1.61717, -1.73788, 5.32027, 8.88357;
  Eigen::MatrixXd sol_UW_exact(2, N);
  sol_UW_exact << -3.29285, -5.76251, 200.931, -1159.4, 3541.58, -0.41613,
      -8.14077, 59.9086, -168.542, -189.822;

  // Compare exact solutions to computed solutions
  ASSERT_LE((sol_LF - sol_LF_exact).norm(), 1e-2);
  ASSERT_LE((sol_UW - sol_UW_exact).norm(), 1e-2);
}

TEST(LinHypDampWaveSys, trackEnergysWavSol) {
  // Some parameters
  double c = 1.;
  double r = 1.;
  double T = 2.;
  unsigned int N = 10;

  // Compute the energies
  auto energies = trackEnergysWavSol(c, r, T, N, UW);

  // Find exact solutions
  Eigen::MatrixXd energies_exact(5, 2);
  energies_exact << 0, 2, .5, 2, 1, 2.75, 1.5, 4, 2, 6;

  // Compare exact to computed solutions
  ASSERT_LE((energies - energies_exact).norm(), 1e-2);
}

}  // namespace LinHypDampWaveSys::test
