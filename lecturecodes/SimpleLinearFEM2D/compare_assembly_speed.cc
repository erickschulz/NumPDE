/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include "SimpleLinearFEM2D.h"
#include "Timer.h"

const double pi = 3.1415926535897;

int main() {
  // configures the number of meshes to try (should be between 1 and 8)
  const int num_meshes = 6;
  // configures the number of repetitions computed per mesh to average over
  const int num_tries = 10;

  // right-hand-side source function f
  FHandle_t f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // Initialize assembly classes
  FESolver solver(f);
  FESolver solver_inefficient(f, 1);
  // 2D array for returning results
  Eigen::Matrix<double, num_meshes, 3> times =
      Eigen::MatrixXd::Zero(num_meshes, 3);
  // Auxiliary object for measuring runtimes
  Timer timer;

  for (int i = 0; i < num_meshes; i++) {
    // load the mesh from file
    std::string mesh_file = "meshes/Square" + std::to_string(i) + ".txt";
    TriaMesh2D mesh(mesh_file);
    // save the number of degrees of freedom
    times(i, 0) = mesh._nodecoords.rows();
    // solve the problem on the mesh num_tries often taking the time for both
    // the efficient and the inefficient method
    Eigen::VectorXd mu;
    for (int j = 0; j < num_tries; j++) {
      std::cout << "Solve " << j + 1 << "/" << num_tries << " of mesh " << i + 1
                << "/" << num_meshes << std::endl;

      // efficient solve
      timer.start();
      mu = solver.Solve(mesh);
      times(i, 1) += timer.elapsed();

      // inefficient solve
      timer.start();
      mu = solver_inefficient.Solve(mesh);
      times(i, 2) += timer.elapsed();
    }
    // average the data and normalize to seconds (measurement is in
    // milliseconds)
    times(i, 1) /= num_tries * 1000.0;
    times(i, 2) /= num_tries * 1000.0;
  }

  return 0;
}
