/**
 * @file symplectictimestepping.cc
 * @brief NPDE homework SymplecticTimestepping
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

constexpr double PI = 3.14159265358979323846;

namespace SymplecticTimestepping {

/* SAM_LISTING_BEGIN_0 */
void sympTimestep(double tau, Eigen::Vector2d& pq_j) {
  // Coefficients of the method
  const Eigen::Vector3d a{2. / 3., -2. / 3., 1.};
  const Eigen::Vector3d b{7. / 24., 3. / 4., -1. / 24.};
  // Single step
#if SOLUTION
  for (int i = 0; i < 3; ++i) {
    pq_j(0) += tau * b(i) * pq_j(1);
    pq_j(1) -= tau * a(i) * pq_j(0);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void sympTimesteppingODETest() {
  int nIter = 7;   // total number of iterations
  unsigned int m;  // number of equidistant steps
  Eigen::Vector2d approx_sol;

  // Initial conditions
  Eigen::Vector2d init_cond;
  init_cond << 0.0, 1.0;  // initial conditions

  // Evaluating the error at the final step between the approx solutions as
  // given by the symplectic method and the exact solution computed from
  // the anlytic formula.
  double errors[nIter];  // errors vector for all approx. sols
  double tau;
  for (int k = 0; k < nIter; k++) {
    m = 10 * std::pow(2, k);
#if SOLUTION
    // Computing approximate solution
    approx_sol = init_cond;
    tau = 2.0 * PI / m;
    for (int i = 0; i < m; i++) {
      sympTimestep(tau, approx_sol);
    }
    // Computing the error in the maximum norm
    errors[k] = std::abs(std::sin(2.0 * PI) - approx_sol[0]) +
                std::abs(std::cos(2.0 * PI) - approx_sol[1]);
#else
    //====================
    // Your code goes here
    //====================
#endif
  }
  // Printing results
  std::cout << "Convergence of Symplectic Time Stepping Method:\n";
  std::cout << "\tM\t\terr(M)\n";
  for (int k = 0; k < nIter; k++) {
    std::cout << "\t" << 10 * std::pow(2, k) << "\t\t" << errors[k]
              << std::endl;
  }
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd simulateHamiltonianDynamics(const Eigen::VectorXd& p0,
                                            const Eigen::VectorXd& q0, double T,
                                            unsigned int M) {
  int n = p0.size();
  Eigen::MatrixXd PQ(2 * n, M + 1);

  // Coefficients of the method
  Eigen::VectorXd a(3);
  a << 2. / 3., -2. / 3., 1.;
  Eigen::VectorXd b(3);
  b << 7. / 24., 3. / 4., -1. / 24.;

  double tau = T / M;
  Eigen::VectorXd pj(p0), qj(q0);
  PQ.col(0) << pj, qj;
#if SOLUTION
  for (int j = 1; j <= M; j++) {   // integrate
    for (int k = 0; k < 3; k++) {  // one step
      // f(q) = - 4 * |q|^2 * q
      pj -= tau * b(k) * 4. * qj.squaredNorm() * qj;
      // g(p) = p
      qj += tau * a(k) * pj;
      PQ.col(j) << pj, qj;
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return PQ;
}
/* SAM_LISTING_END_3 */

}  // namespace SymplecticTimestepping
