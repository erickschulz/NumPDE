/**
 * @ file stokespipeflow_main.cc
 * @ brief NPDE homework StokesPipeFlow MAIN FILE
 * @ author Wouter Tonnon
 * @ date May 2024
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/base/lehrfem_info.h>
#include <lf/base/ref_el.h>

#include <iostream>
#include <memory>
#include <vector>

#include "stokespipeflow.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "NumPDE homework problem StokesPipeFlow by Wouter Tonnon (R)\n";
  lf::base::LehrFemInfo::PrintInfo(std::cout);

  // Test of convergence on the unit square with manufactured solution
  //  StokesPipeFlow::testCvgTaylorHood(4);
  // Compute on "realistic" geometry
  //  std::cout << "Writing vtk-File" << std::endl;
  StokesPipeFlow::visualizeTHPipeFlow("meshes/pipe3.msh", "out.vtk");
  std::vector<std::pair<double, double>> p_diss{};
  double p_diss_vol;
  double p_diss_bd;
  for (char no : {'1', '2', '3'}) {
    std::string meshfile = std::string("meshes/pipe") + no + ".msh";
    std::cout << "Reading mesh from file " << meshfile << std::endl;
    p_diss_vol = StokesPipeFlow::computeDissipatedPower(meshfile.c_str());
    p_diss_bd = StokesPipeFlow::computeDissipatedPowerBd(meshfile.c_str());
    std::cout << "Dissipated power: volume formula = " << p_diss_vol
              << ", boundary formula = " << p_diss_bd << std::endl;
    p_diss.emplace_back(p_diss_vol, p_diss_bd);
  }
  std::cout << "level "
            << " p_vol "
            << " p_bd "
            << " D(p_vol) "
            << " D(p_bd)\n";
  for (int k = 0; k < p_diss.size(); ++k) {
    std::cout << "l = " << k << ": p_vol = " << p_diss[k].first
              << ", p_bd = " << p_diss[k].second
              << ", D(p_vol) = " << p_diss[k].first - p_diss_vol
              << ", D(p_bd) = " << p_diss[k].second - p_diss_vol << std::endl;
  }

  return 0;
}
