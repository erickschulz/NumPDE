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

  /* Your code goes here */

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

  /* Your code goes here */

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

  /* Your code goes here */
  Eigen::MatrixXd sol = Eigen::MatrixXd::Zero(2, N);
  // Return
  return sol;
}
/* SAM_LISTING_END_4 */

Eigen::Matrix<double, Eigen::Dynamic, 2> trackEnergysWavSol(
    double c, double r, double T, unsigned int N, Fluxtype nf_selector);

}  // namespace LinHypDampWaveSys

#endif
