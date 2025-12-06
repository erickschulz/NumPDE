/**
 * @ file stokesstabp1fem_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <iostream>

#include "stokesstabp1fem.h"

int main(int argc, char** argv) {
  std::cout << "NumPDE homework problem StokesStabP1FEM\n";
  std::cout << "Created by R. Hiptmair, May 2025\n";
  lf::base::LehrFemInfo::PrintInfo(std::cout);

  unsigned int refsteps = 4;
  if (argc > 1) {
    refsteps = atoi(argv[1]);
  }
  if (refsteps > 0) {
    StokesStabP1FEM::testCvgP1StabFEM(6);
  }
  return 0;
}
