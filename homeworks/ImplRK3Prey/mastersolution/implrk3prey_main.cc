/**
 * @file implrk3prey_main.cc
 * @brief NPDE homework ImplRK3Prey code
 * @author Unknown, Oliver Rietmann
 * @date 29.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "implrk3prey.h"

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Definition of coefficients in Butcher scheme
  unsigned int s = 2;
  Eigen::MatrixXd A(s, s);
  Eigen::VectorXd b(s);
  // What method is this?
  A << 5. / 12., -1. / 12., 3. / 4., 1. / 4.;
  b << 3. / 4., 1. / 4.;
  // Initialize implicit RK with Butcher scheme
  ImplRK3Prey::implicitRKIntegrator integrator(A, b);
  /* SAM_LISTING_END_0 */

  /* SAM_LISTING_BEGIN_1 */
  // Coefficients and handle for prey/predator model
  double alpha1 = 3.;
  double alpha2 = 2.;
  double beta1 = 0.1;
  double beta2 = 0.1;
  auto f = [&alpha1, &alpha2, &beta1, &beta2](const Eigen::VectorXd& y) {
    auto temp = y;
    temp(0) *= alpha1 - beta1 * y(1);
    temp(1) *= -alpha2 + beta2 * y(0);
    return temp;
  };

  auto Jf = [&alpha1, &alpha2, &beta1, &beta2](const Eigen::VectorXd& y) {
    Eigen::Matrix2d temp;
    temp << alpha1 - beta1 * y(1), -beta1 * y(0), beta2 * y(1),
        -alpha2 + beta2 * y(0);
    return temp;
  };

  // Final time for model
  double T = 10.;

  // Initial value for model
  Eigen::Vector2d y0(100.0, 5.0);

  // Array of number of steps (for convergence study)
  std::vector<unsigned int> N = {128,  256,  512,   1024,  2048,
                                 4096, 8192, 16384, 32768, 65536};

  // Exact value y(10) at final time T = 10 (approximated)
  Eigen::Vector2d yex(0.319465882659820, 9.730809352326228);
  // Start convergence study
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
  double errold;
  for (unsigned int i = 0; i < N.size(); ++i) {
    auto res = integrator.solve(f, Jf, T, y0, N[i]);
    double err = (res.back() - yex).norm();
    std::cout << std::setw(15) << N[i] << std::setw(15) << err;
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(errold / err);
    }
    errold = err;
    std::cout << std::endl;
  }
  /* SAM_LISTING_END_1 */

  return 0;
}
