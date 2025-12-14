/**
 * @ file magdiffwire_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Ralf Hiptmair
 * @ date June 2025
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <iostream>

#include "magdiffwire.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "NumPDE homework problem MagDiffWire\n";
  std::cout << "Coded by R. Hiptmair, June 2025\n";
  lf::base::LehrFemInfo::PrintInfo(std::cout);

  MagDiffWire::testCvgMagDiffWire(4);

  return 0;
}

/*
Original version of sdirkMagDiffWire
N                   L2 err         H1s err
---------------------------------------------
40        0.0469          0.46
139       0.0524          0.502
517       0.0531          0.508
1993      0.0532          0.508
7825      0.0532          0.508

Alternative version:

N                   L2 err         H1s err
---------------------------------------------
40        0.517           2.32
139       0.523           2.32
517       0.523           2.33
1993      0.523           2.33
7825      0.523           2.33

Matrix version



 */
