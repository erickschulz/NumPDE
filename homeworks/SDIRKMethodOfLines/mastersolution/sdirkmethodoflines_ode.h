/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <iostream>
#include <vector>

namespace SDIRKMethodOfLines {

std::vector<double> sdirk2SteppingLinScalODE(unsigned int m);

void sdirk2ScalarODECvTest();

}  // namespace SDIRKMethodOfLines
