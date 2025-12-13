/**
 * @ file leastsquaresadvection_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Ralf Hiptmair
 * @ date July 2024
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/uscalfe/fe_space_lagrange_o1.h>

#include <cmath>
#include <iostream>

#include "leastsquaresadvection.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout
      << "LehrFEM++-based NumPDE homework problem LeastSquaresAdvection\n\n";

  {
    std::cout << "\n >>> Setting (I): Convergence for non-smmoth solution\n";
    auto g = [](Eigen::Vector2d x) -> double {
      LF_ASSERT_MSG(
          (x[0] >= 0.0) and (x[0] <= 1.0) and (x[1] >= 0.0) and (x[1] <= 1.0),
          "x must lie inside the unit sqaure!");
      if (x[1] < 1E-8) {
        const double s = std::sin(M_PI * x[0]);
        return s * s;
      }
      return 0.0;
    };
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO1<double>>(g, 1.0, 3);
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO2<double>>(g, 1.0, 3);
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO3<double>>(g, 1.0, 3);
  }

  {
    std::cout << " >>> Setting (II): Convergence for smooth solution\n";
    auto g = [](Eigen::Vector2d x) -> double {
      LF_ASSERT_MSG(
          (x[0] >= 0.0) and (x[0] <= 1.0) and (x[1] >= 0.0) and (x[1] <= 1.0),
          "x must lie inside the unit sqaure!");
      const double xi = Eigen::Vector2d(-1.0, 2.0).dot(x);
      return std::cos(M_PI * xi);
    };
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO1<double>>(g, 0.0, 3);
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO2<double>>(g, 0.0, 3);
    LeastSquaresAdvection::testCVGLSQAdvectionReaction<
        lf::uscalfe::FeSpaceLagrangeO3<double>>(g, 0.0, 3);
  }

  return 0;
}
