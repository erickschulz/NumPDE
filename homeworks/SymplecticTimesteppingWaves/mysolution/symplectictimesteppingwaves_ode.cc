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
  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
}

/* SAM_LISTING_END_2 */

}  // namespace SymplecticTimesteppingWaves
