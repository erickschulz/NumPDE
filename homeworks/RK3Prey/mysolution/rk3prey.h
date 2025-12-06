/**
 * @file rk3prey_main.cc
 * @brief NPDE homework RK3Prey code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <vector>

namespace RK3Prey {

// Butcher Tableau based Runge-Kutta explicit solver for autonomous ODEs
/* SAM_LISTING_BEGIN_0 */
class RKIntegrator {
 public:
  RKIntegrator(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    //====================
    // Your code goes here
    //====================
  }

  // Explicit Runge-Kutta numerical integrator
  template <class Function>
  std::vector<Eigen::VectorXd> solve(Function&& f, double T,
                                     const Eigen::VectorXd& y0, int M) const;

 private:
  //====================
  // Your code goes here
  //====================
};
/* SAM_LISTING_END_0 */

/* Solves an autonomous ODE y' = f(y), y(0) = y0, using a
 * RK scheme from the Butcher tableau provided by the
 * constructor. Performs N equidistant steps up to time T */
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
std::vector<Eigen::VectorXd> RKIntegrator::solve(Function&& f, double T,
                                                 const Eigen::VectorXd& y0,
                                                 int M) const {
  int dim = y0.size();  // dimension
  double h = T / M;     // step size
  std::vector<Eigen::VectorXd> sol;
  sol.reserve(M + 1);

  //====================
  // Your code goes here
  //====================
  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace RK3Prey
