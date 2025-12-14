#include "ordnotall.h"

namespace OrdNotAll {

/* SAM_LISTING_BEGIN_1 */
template <class Function>
void testCvgRKSSM(const Function& f, double T, double y0,
                  const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
  // Helper object carrying out the actual explicit RK-SSM
  RKIntegrator<double> rk(A, b);
  // Vector for collecting errors
  std::vector<double> error(15);
  // TODO: output error and order of the method
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void cmpCvgRKSSM() {
  // Construct data for Butcher schemes
  Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(1, 1);
  Eigen::VectorXd b1(1);
  b1 << 1;

  Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(2, 2);
  A2(1, 0) = 1;
  Eigen::VectorXd b2(2);
  b2 << .5, .5;

  Eigen::MatrixXd A3 = Eigen::MatrixXd::Zero(3, 3);
  A3(1, 0) = .5;
  A3(2, 0) = -1;
  A3(2, 1) = 2;
  Eigen::VectorXd b3(3);
  b3 << 1. / 6, 2. / 3, 1. / 6;

  Eigen::MatrixXd A4 = Eigen::MatrixXd::Zero(4, 4);
  A4(1, 0) = .5;
  A4(2, 1) = .5;
  A4(3, 2) = 1;
  Eigen::VectorXd b4(4);
  b4 << 1. / 6, 1. / 3, 1. / 3, 1. / 6;

  // TODO: call testCvgRKSSM for all combinations of ODE and RK methods
}
/* SAM_LISTING_END_2 */

}  // namespace OrdNotAll
