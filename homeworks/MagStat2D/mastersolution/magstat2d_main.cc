/**
 * @ file magstat2d_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Ralf Hiptmair
 * @ date June 2025
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <iostream>

#include "magstat2d.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "NumPDE homework problem MagStat2D\n";
  std::cout << "Created by R. Hiptmair, June 2025\n";
  lf::base::LehrFemInfo::PrintInfo(std::cout);

  // Empiric convergence study
  MagStat2D::testCvgMagStat2DWhitneyFEM(4);

  return 0;
}
