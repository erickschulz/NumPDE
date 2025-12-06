/**
 * @ file transformedconslaw_main.cc
 * @ brief NPDE exam problem
 * @ author Oliver Rietmann
 * @ date 05.08.2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

#include "systemcall.h"
#include "transformedconslaw.h"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  // Setup
  TRFCL::NonStdCauchyProblemCL prb;
  int M = 500;
  int N = 500;

  // Compute cell points and inital data for plot
  std::pair<double, double> limits = prb.domain();
  Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, limits.first, limits.second);
  auto z0 = [&prb](double y) { return prb.z0(y); };
  Eigen::VectorXd zeta0 = x.unaryExpr(z0);

  // Print conserved quantity
  auto rho = [&prb](double z) { return prb.rho(z); };
  auto rec = [rho, limits](double t, const Eigen::VectorXd& zeta) -> void {
    Eigen::VectorXd weighted_samples = zeta.unaryExpr(rho) / (zeta.size() - 1);
    double rho_z_integral =
        weighted_samples.sum() * (limits.second - limits.first);
    std::cout << "Integral of rho(z(x,t)) at time t = " << t << ":\t"
              << rho_z_integral << std::endl;
  };

  // Compute solution at final time
  std::cout << std::fixed;
  Eigen::VectorXd zetaT = solveCauchyPrb(M, N, prb, rec);

  // Write x, zeta0 and zetaT to .csv file
  std::ofstream solution_file;
  solution_file.open("solution.csv");
  solution_file << x.transpose().format(CSVFormat) << std::endl;
  solution_file << zeta0.transpose().format(CSVFormat) << std::endl;
  solution_file << zetaT.transpose().format(CSVFormat) << std::endl;
  solution_file.close();
  std::cout << "Generated solution.csv" << std::endl;

  // Plot from .csv file using plot.py
  systemcall::execute("python3 scripts/plot.py solution.csv solution.eps");

  return 0;
}
