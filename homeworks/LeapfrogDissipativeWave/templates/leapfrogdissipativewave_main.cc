/**
 * @ file leapfrogdissipativewave_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author R. Hiptmair
 * @ date July 2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <boost/program_options.hpp>
#include <iostream>

#include "leapfrogdissipativewave.h"

int main(int argc, char** argv) {
  std::cout << "Mastersolution problem LeapfrogDissipativeWave" << std::endl;
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "Help message")
  ("scalar,s", "Scalar test case")  
  ("levels,l", po::value<unsigned int>()->default_value(5), "no of refinement levels")
  ("M0,m", po::value<unsigned int>()->default_value(100), "no timesteps of coarsest level")
  ("Mfactor,f", po::value<unsigned int>()->default_value(2), "increase factor for M");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
  } else {
    if (vm.count("scalar") > 0) {
      LeapfrogDissipativeWave::testDissipativeLeapfrog();
      LeapfrogDissipativeWave::convergenceScalarTestProblem();
    } else {
      unsigned int reflevel = vm["levels"].as<unsigned int>();
      unsigned int M0 = vm["M0"].as<unsigned int>();
      unsigned int Mfac = vm["Mfactor"].as<unsigned int>();
      LeapfrogDissipativeWave::convergenceDissipativeLeapfrog(reflevel, 1.0, M0,
                                                              Mfac);
    }
  }
  return 0;
}
