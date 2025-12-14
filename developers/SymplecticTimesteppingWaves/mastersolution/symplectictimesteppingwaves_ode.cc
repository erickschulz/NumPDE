/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves_ode.h"

namespace SymplecticTimesteppingWaves {

/* SAM_LISTING_BEGIN_1 */
void sympTimestep(double tau, Eigen::Vector2d& pq_j) {
#if SOLUTION
  // Coefficients of the method
  Eigen::VectorXd a(3);
  a << 2. / 3., -2. / 3., 1.;
  Eigen::VectorXd b(3);
  b << 7. / 24., 3. / 4., -1. / 24;

  // one step method
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
/* SAM_LISTING_END_1 */

Eigen::Vector2d sympTimesteppingHarmonicOscillatorODE(unsigned int m) {
  Eigen::Vector2d approx_sol;
  approx_sol << 0.0, 1.0;  // initial conditions
  double tau = 2.0 * M_PI / m;
  for (int i = 0; i < m; i++) {
    sympTimestep(tau, approx_sol);
  }
  return approx_sol;
}

/* SAM_LISTING_BEGIN_2 */
void sympTimesteppingODETest() {
#if SOLUTION
  int nIter = 9;   // total number of iterations
  unsigned int m;  // number of equidistant steps

  // Evaluating the error at the final step between the approx solutions as
  // given by the symplectic method and the exact solution computed from
  // the anlytic formula
  Eigen::Vector2d approx_sol;
  double errors[nIter];  // errors vector for all approx. sols
  for (int k = 0; k < nIter; k++) {
    m = 10 * std::pow(2, k);
    // Computing approximate solution
    approx_sol = sympTimesteppingHarmonicOscillatorODE(m);
    // Computing the error in the maximum norm
    errors[k] = std::abs(sin(2.0 * M_PI) - approx_sol[0]) +
                std::abs(cos(2.0 * M_PI) - approx_sol[1]);
  }
  // Computing rates of convergence
  double rates[nIter - 1];
  double avg_rate = 0.0;

  for (int k = 0; k < nIter - 1; k++) {
    rates[k] = log2(errors[k] / errors[k + 1]);
    avg_rate += rates[k];
  }
  avg_rate = avg_rate / (nIter - 1);

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "     Convergence of Symplectic Time Stepping Method      "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| Nsteps"
            << "\t| error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << k << "\t"
              << "\t|" << 10 * std::pow(2, k) << "\t\t|" << errors[k];
    if (k > 0) {
      std::cout << "\t|" << rates[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "Average rate of convergence: " << avg_rate << "\n" << std::endl;
#else
  //====================
  // Your code goes here
  //====================
#endif
}

/* SAM_LISTING_END_2 */

}  // namespace SymplecticTimesteppingWaves
