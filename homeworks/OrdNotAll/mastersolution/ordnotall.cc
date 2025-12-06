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
  // Vector for storing rate estimates
  std::vector<double> order(14);

  double sum = 0;
  int count = 0;
  bool test = true;
  // Reference numerical solution obtained with $2^{15}$ timesteps
  std::vector<double> y_exact = rk.solve(f, T, y0, std::pow(2, 15));

  for (int k = 0; k < 12; k++) {
    // Number of timesteps
    int M = std::pow(2, k + 1);
    // Solve IVP
    std::vector<double> y1 = rk.solve(f, T, y0, M);
    // Error at final time
    error[k] = std::abs(y1[M] - y_exact[std::pow(2, 15)]);

    std::cout << std::left << std::setfill(' ') << std::setw(3)
              << "M = " << std::setw(7) << M << std::setw(8)
              << "Error = " << std::setw(13) << error[k];

    if (error[k] < 1e-14) {
      test = false;
    }
    if (k > 0 && test) {
      order[k - 1] = std::log2(error[k - 1] / error[k]);
      std::cout << std::left << std::setfill(' ') << std::setw(10)
                << "Approximated order = " << order[k - 1] << std::endl;
      sum += order[k - 1];
      count = k;
    } else
      std::cout << std::endl;
  }
  std::cout << "Average approximated order = " << sum / count << std::endl
            << std::endl;
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

  // First ODE
  std::cout << std::endl
            << "1. ODE y' = (1-y)y, y(0)=.5" << std::endl
            << std::endl;
  double T = 1;
  auto f = [](double y) { return (1. - y) * y; };
  double y0 = .5;

  std::cout << "Explicit Euler" << std::endl << std::endl;
  testCvgRKSSM(f, T, y0, A1, b1);
  std::cout << "Trapezoidal rule" << std::endl << std::endl;
  testCvgRKSSM(f, T, y0, A2, b2);
  std::cout << "RK order 3" << std::endl << std::endl;
  testCvgRKSSM(f, T, y0, A3, b3);
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  testCvgRKSSM(f, T, y0, A4, b4);

  // Second ODE
  std::cout << std::endl
            << "2. ODE y' = |1.1 - y| + 1, y(0)=1" << std::endl
            << std::endl;
  auto f2 = [](double y) { return std::abs(1.1 - y) + 1.; };
  y0 = 1;

  std::cout << "Explicit Euler" << std::endl << std::endl;
  testCvgRKSSM(f2, T, y0, A1, b1);
  std::cout << "Trapezoidal rule" << std::endl << std::endl;
  testCvgRKSSM(f2, T, y0, A2, b2);
  std::cout << "RK order 3" << std::endl << std::endl;
  testCvgRKSSM(f2, T, y0, A3, b3);
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  testCvgRKSSM(f2, T, y0, A4, b4);
}
/* SAM_LISTING_END_2 */

}  // namespace OrdNotAll
