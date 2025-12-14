/**
 * @file finitevolumesineconslaw_main.cc
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <fstream>
#include <iostream>

#include "finitevolumesineconslaw.h"
#include "systemcall.h"

using namespace FiniteVolumeSineConsLaw;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  /* SAM_LISTING_BEGIN_1 */
  unsigned int N = 600;
  unsigned int M = 200;
  double h = 12.0 / N;
  Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, -6.0 + 0.5 * h, 6.0 - 0.5 * h);

  std::ofstream file;

  // without reaction term
  Eigen::VectorXd ufinal = solveSineConsLaw(&sineClawRhs, N, M);
#if SOLUTION
  file.open("ufinal.csv");
  file << x.transpose().format(CSVFormat) << std::endl;
  file << ufinal.transpose().format(CSVFormat) << std::endl;
  file.close();

  std::cout << "Generated ufinal.csv" << std::endl;
  systemcall::execute("python3 scripts/plot.py ufinal.csv ufinal.eps");
#else
  //====================
  // Your code goes here
  // Use std::ofstream to write the solution to
  // the file "ufinal.csv". To plot this
  // file you may uncomment the following line:
  // systemcall::execute("python3 scripts/plot.py ufinal.csv
  // ufinal.eps");
  //====================
#endif
  /* SAM_LISTING_END_1 */

  // with reaction term: -c * u(x, t), where c = 1.0
  double c = 1.0;
  auto bind_c = [c](const Eigen::VectorXd& mu) {
    return sineClawReactionRhs(mu, c);
  };
  Eigen::VectorXd ufinal_reaction = solveSineConsLaw(bind_c, N, M);
#if SOLUTION
  file.open("ufinal_reaction.csv");
  file << x.transpose().format(CSVFormat) << std::endl;
  file << ufinal_reaction.transpose().format(CSVFormat) << std::endl;
  file.close();

  std::cout << "Generated ufinal_reaction.csv" << std::endl;
  systemcall::execute(
      "python3 scripts/plot.py ufinal_reaction.csv "
      "ufinal_reaction.eps");
#else
  //====================
  // Your code goes here
  // Use std::ofstream to write the solution to
  // the file "ufinal_reaction.csv". To plot this
  // file you may uncomment the following line:
  // systemcall::execute("python3 scripts/plot.py ufinal_reaction.csv
  // ufinal_reaction.eps");
  //====================
#endif

  // Finding the optimal timestep (no reaction term)
  unsigned int M_small = findTimesteps();
  std::cout
      << "The smallest number of timesteps keeping the solution in [0, 2] is "
      << M_small << std::endl;

  return 0;
}
