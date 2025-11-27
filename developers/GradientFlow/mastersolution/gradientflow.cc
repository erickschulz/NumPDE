/**
 * @file gradientflow.h
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "gradientflow.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace GradientFlow {

/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd ButcherMatrix() {
  Eigen::MatrixXd A(6, 5);
  // clang-format off
  A <<       0.25,          0.,       0.,       0.,   0.,
    0.5,        0.25,       0.,       0.,   0.,
    17./50.,     -1./25.,     0.25,       0.,   0.,
    371./1360., -137./2720., 15./544.,     0.25,   0.,
    25./24.,    -49./48., 125./16., -85./12., 0.25,
    25./24.,    -49./48., 125./16., -85./12., 0.25;
  // clang-format on
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::VectorXd> SolveGradientFlow(const Eigen::VectorXd& d,
                                               double lambda,
                                               const Eigen::VectorXd& y0,
                                               double T, unsigned int M) {
  // initialize solution vector
  std::vector<Eigen::VectorXd> sol(M + 1, Eigen::VectorXd::Zero(y0.size()));

#if SOLUTION
  // Define the right hand side of the ODE y' = f(y), and the Jacobian of f.
  auto f = [d, lambda](const Eigen::VectorXd& y) {
    Eigen::VectorXd val =
        -2. * std::cos(y.squaredNorm()) * y - 2. * lambda * d.dot(y) * d;
    return val;
  };
  auto df = [d, lambda](const Eigen::VectorXd& y) {
    int dim = y.size();
    Eigen::MatrixXd term1 = 4. * std::sin(y.squaredNorm()) * y * y.transpose();
    Eigen::MatrixXd term2 =
        -2. * std::cos(y.squaredNorm()) * Eigen::MatrixXd::Identity(dim, dim);
    Eigen::MatrixXd term3 = -2. * lambda * d * d.transpose();
    Eigen::MatrixXd dfval = term1 + term2 + term3;
    return dfval;
  };

  // Split the interval [0,T] into M intervals of size h.
  double h = T / M;
  Eigen::VectorXd y = y0;
  sol[0] = y;
  // Evolve up to time T:
  for (int i = 1; i <= M; i++) {
    y = DiscEvolSDIRK(f, df, y, h);
    sol[i] = y;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace GradientFlow
