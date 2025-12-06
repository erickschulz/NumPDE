/**
 * @ file LinearFE1D_main.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher, Erick Schulz
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

/* This problem does not rely on Lehrfempp. It is solely based on the
library Eigen. As the computational domain is simply the 1D interval
(0,1), the mesh is stored as a vector whose entries are the
coordinates of the nodes (i.e. the distance of the node from origin 0.0).
*/

#include "linearfe1d.h"

int main() {
  // Creating a 1D mesh of the interval (0,1)
  int N = 100;                  // nb. of cells
  Eigen::VectorXd mesh(N + 1);  // nb. of nodes
  // Nodes are equally spaced
  for (int i = 0; i < N + 1; i++) {
    mesh[i] = i * (1.0 / N);
  }

  // Constant and variable parameters
  auto identity = [](double x) { return x; };
  auto const_one = [](double x) { return 1.0; };

  // Solving the BVPs
  Eigen::VectorXd uA = LinearFE1D::solveA(mesh, identity, const_one);
  Eigen::VectorXd uB = LinearFE1D::solveB(mesh, identity, const_one, 0.1, 0.5);
  Eigen::VectorXd uC = LinearFE1D::solveC(mesh, identity, identity);

  // PRINTING RESULTS TO.csv FILE
  // Defining CSV output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream ofs;
  // Printing results to file for problem (A)
  std::string filename = "uA.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uA.format(CSVFormat);
    std::cout << "Components of uA were written to uA.csv" << std::endl;
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uA.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (B)
  filename = "uB.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uB.format(CSVFormat);
    std::cout << "Components of uB were written to uB.csv" << std::endl;
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uB.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (C)
  filename = "uC.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uC.format(CSVFormat);
    std::cout << "Components of uC were written to uC.csv" << std::endl;
  }
  ofs.close();
}
