/**
 * @ file boundarylength_main.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <iostream>
#include <utility>

#include "boundarylength.h"

using namespace LengthOfBoundary;

/* SAM_LISTING_BEGIN_1 */
int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string file_name(argv[1]);
    std::pair<double, double> result = measureDomain(file_name);

    std::cout << "The area of the domain is: " << result.first << std::endl;
    std::cout << "The length of the boundary is: " << result.second
              << std::endl;
  } else {
    std::cout << "Correct usage: Call this executable and supply the name of a "
                 "msh-file as a command line argument "
              << std::endl;
  }
}
/* SAM_LISTING_END_1 */
