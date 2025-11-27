/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */
#include <iostream>

#include "SimpleLinearFEM2D.h"

// define pi
const double pi = 3.1415926535897;

int main() {
  // load mesh
  TriaMesh2D square_mesh("meshes/Square5.txt");
  std::cout << "Mesh loaded " << std::endl;
  std::cout << "Mesh info: " << square_mesh._nodecoords.rows() << " vertices,  "
            << square_mesh._elements.rows() << " elements" << std::endl;

  // source function f
  FHandle_t f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // solve the system
  FESolver solver(f);
  Eigen::VectorXd solution = solver.Solve(square_mesh);

  return 0;
}
