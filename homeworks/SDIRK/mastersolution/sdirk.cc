/**
 * @file sdirk.cc
 * @brief NPDE homework SDIRK code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "sdirk.h"

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.h"

namespace SDIRK {

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d SdirkStep(const Eigen::Vector2d& z0, double h, double gamma) {
  Eigen::Vector2d res;
  // Compute one timestep of the SDIRK implicit RK-SSM for the linear ODE
  // Matrix A for evaluation of f
  Eigen::Matrix2d A;
  A << 0., 1., -1., -1.;
  // Precompute and reuse factorization
  auto A_lu = (Eigen::Matrix2d::Identity() - h * gamma * A).partialPivLu();
  Eigen::Vector2d az = A * z0;

  // Increments according to \prbeqref{eq:ies}
  Eigen::Vector2d k1 = A_lu.solve(az);
  Eigen::Vector2d k2 = A_lu.solve(az + h * (1 - 2 * gamma) * A * k1);

  // Updated state
  res = z0 + h * 0.5 * (k1 + k2);
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> SdirkSolve(const Eigen::Vector2d& z0,
                                        unsigned int M, double T,
                                        double gamma) {
  // Solution vector
  std::vector<Eigen::Vector2d> res(M + 1);
  // Solve the ODE with uniform timesteps using the SDIRK method
  // Equidistant step size
  const double h = T / M;
  // Push initial data
  res[0] = z0;
  // Main loop
  for (unsigned int i = 1; i <= M; ++i) {
    res[i] = SdirkStep(res[i - 1], h, gamma);
  }
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double CvgSDIRK() {
  double conv_rate;
  // Study the convergence rate of the method.
  // Initial data z0 = [y(0), y'(0)]
  Eigen::Vector2d z0;
  z0 << 1, 0;
  // Final time
  const double T = 10;
  // Parameter
  const double gamma = (3. + std::sqrt(3.)) / 6.;
  // Mesh sizes
  Eigen::ArrayXd err(10);
  Eigen::ArrayXd M(10);
  M << 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240;

  // Exact solution (only y(t)) given z0 = [y(0), y'(0)] and t
  auto yex = [&z0](double t) {
    return 1. / 3. * std::exp(-t / 2.) *
           (3. * z0(0) * std::cos(std::sqrt(3.) * t / 2.) +
            std::sqrt(3.) * z0(0) * std::sin(std::sqrt(3.) * t / 2.) +
            2. * std::sqrt(3.) * z0(1) * std::sin(std::sqrt(3.) * t / 2.));
  };

  // Store old error for rate computation
  double errold = 0;
  std::cout << std::setw(15) << "m" << std::setw(15) << "maxerr"
            << std::setw(15) << "rate" << std::endl;
  // Loop over all meshes
  for (unsigned int i = 0; i < M.size(); ++i) {
    int m = M(i);
    // Get solution
    auto sol = SdirkSolve(z0, m, T, gamma);
    // Compute error
    err(i) = std::abs(sol.back()(0) - yex(T));

    // Print table
    std::cout << std::setw(15) << m << std::setw(15) << err(i);
    if (i > 0) std::cout << std::setw(15) << std::log2(errold / err(i));
    std::cout << std::endl;

    // Store old error
    errold = err(i);
  }

  Eigen::VectorXd coeffs = polyfit(M.log(), err.log(), 1);
  conv_rate = -coeffs(0);
  return conv_rate;
}
/* SAM_LISTING_END_2 */

}  // namespace SDIRK
