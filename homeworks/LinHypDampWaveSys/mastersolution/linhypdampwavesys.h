/**
 * @file linhypdampwavesys.h
 * @brief NPDE homework LinHypDampWaveSys code
 * @author Wouter Tonnon
 * @date May 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef LinHypDampWaveSys_H_
#define LinHypDampWaveSys_H_

#include <Eigen/Core>
#include <iomanip>
#include <iostream>

namespace LinHypDampWaveSys {
enum Fluxtype { LF, UW };

void visWavSol(double c, double r, double T, unsigned int N,
               Fluxtype nf_selector);

/* SAM_LISTING_BEGIN_1 */
template <typename FFunctor, typename gFunctor>
Eigen::MatrixXd linsysoderhs(const Eigen::MatrixXd& mu, double h,
                             FFunctor&& numflux, gFunctor&& g) {
  unsigned n = mu.cols();                             // length of state vector
  Eigen::MatrixXd fd = Eigen::MatrixXd::Zero(2, n);   // flux difference vector
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2, n);  // return vector

  // constant continuation of data
  fd.col(0) = numflux(mu.col(0), mu.col(1)) - numflux(mu.col(0), mu.col(0));
  // Fluxes on the interior
  for (unsigned j = 1; j < n - 1; ++j) {
    fd.col(j) = numflux(mu.col(j), mu.col(j + 1)) -
                numflux(mu.col(j - 1), mu.col(j));  // see \eqref{eq:2pcf}
  }
  // constant continuation of data!
  fd.col(n - 1) = numflux(mu.col(n - 1), mu.col(n - 1)) -
                  numflux(mu.col(n - 2), mu.col(n - 1));

  // Compute the rhs
  for (unsigned j = 0; j < n; ++j) {
    rhs.col(j) = -fd.col(j) / h + g(mu.col(j));
  }

  // Efficient thanks to return value optimization (RVO)
  return rhs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
template <typename FFunctor, typename gFunctor, typename u0Functor,
          typename RECORDER = std::function<void(const Eigen::MatrixXd&)>>
Eigen::MatrixXd fvEulLinSys(
    double a, double b, double T, unsigned int N, unsigned int M,
    u0Functor&& u0, FFunctor& numflux, gFunctor&& f,
    RECORDER&& recorder = [](const Eigen::MatrixXd&) -> void {}) {
  // Determine constants of scheme
  double dt = T / M;
  double h = (b - a) / N;
  Eigen::MatrixXd mu(2, N);  // return vector

  // Set initial conditions
  for (unsigned j = 0; j < N - 1; ++j) {
    mu.col(j) = u0(a + h / 2. + j * h);
  }

  // Solve the ODE
  recorder(mu);
  for (double t = 0; t < T - dt / 2.; t += dt) {
    Eigen::MatrixXd rhs = linsysoderhs(mu, h, numflux, f);
    mu += dt * rhs;
    recorder(mu);
  }

  // Return
  return mu;
}

/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename u0Functor,
          typename RECORDER = std::function<void(const Eigen::MatrixXd&)>>
Eigen::Matrix<double, 2, Eigen::Dynamic> ev1DampWave(
    double c, double r, double T, unsigned int N, u0Functor&& u0,
    Fluxtype nf_selector,
    RECORDER&& recorder = [](const Eigen::MatrixXd&) -> void {}) {
  // Some parameters for the numerical scheme
  double maxspeed = c;
  double a = -T * maxspeed;
  double b = 1 + T * maxspeed;
  double h = (b - a) / N;
  double dt = h / maxspeed;  //*(b - a) / ((N * maxspeed));
  unsigned int M = std::ceil(T / dt);

  // Define numerical flux
  auto numflux = [c, nf_selector](Eigen::VectorXd mu1,
                                  Eigen::VectorXd mu2) -> Eigen::VectorXd {
    Eigen::Vector2d out;
    if (nf_selector == UW) {  // Upwind
      out(0) = .5 * (c * mu1(0) - c * c * mu1(1)) +
               .5 * (-c * mu2(0) - c * c * mu2(1));
      out(1) = .5 * (-1 * mu1(0) + c * mu1(1)) + .5 * (-mu2(0) - c * mu2(1));
    } else {  // Lax-Friedrichs
      out(0) = -.5 * c * c * (mu1(1) + mu2(1)) - .5 * c * (mu2(0) - mu1(0));
      out(1) = -.5 * (mu1(0) + mu2(0)) - .5 * c * (mu2(1) - mu1(1));
    }
    return out;
  };

  // The RHS
  auto g = [r](Eigen::VectorXd u) -> Eigen::VectorXd {
    Eigen::Vector2d out;
    out(0) = -r * u(0);
    out(1) = 0.;
    return out;
  };

  // Solve the system
  Eigen::MatrixXd sol = fvEulLinSys(-T * maxspeed, 1 + T * maxspeed, T, N, M,
                                    u0, numflux, g, recorder);
  // Return
  return sol;
}
/* SAM_LISTING_END_4 */

Eigen::Matrix<double, Eigen::Dynamic, 2> trackEnergysWavSol(
    double c, double r, double T, unsigned int N, Fluxtype nf_selector);

}  // namespace LinHypDampWaveSys

#endif
