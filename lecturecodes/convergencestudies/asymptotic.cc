/**
 * @file asymptotic.cc
 * @brief Creates convergence plots for experiment 3.2.3.12
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES  // NOLINT(bugprone-reserved-identifier)

#include <lf/quad/gauss_quadrature.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace po = boost::program_options;

// Code for empiric exploration of asymptotic convergence of norms of the
// discretization error for 1D finite element discretization of a 2-point BVP.
// Special example with a highly oscillatory solution
int main(int argc, char* argv[]) {
  std::cout << "Run with --help for usage information" << std::endl;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "Display this help message")
  ("output,o", po::value<std::string>()->default_value("output.csv"), "Name of the output file")
  ("M_max,M", po::value<int>()->default_value(500), "Maximum number of cells")
  ("dM,m", po::value<int>()->default_value(5), "Increment in M")
  ("num_quad_points,n", po::value<int>()->default_value(2), "Number of points for numerical quadrature");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }
  const int M_max = vm["M_max"].as<int>();
  const int dM = vm["dM"].as<int>();
  const int num_quad_points = vm["num_quad_points"].as<int>();
  const std::string output_file = vm["output"].as<std::string>();
  const auto [quad_points, quad_weights] = lf::quad::GaussLegendre(2);

  // Load function
  const auto f = [](double x) {
    return 10000 * M_PI * M_PI * x * x * std::sin(50 * M_PI * x * x) -
           100 * M_PI * std::cos(50 * M_PI * x * x);
  };
  // Analytic solution (highly oscillatory)
  const auto u = [](double x) { return std::sin(50 * M_PI * x * x); };
  // Gradient of analytic solution
  const auto u_grad = [](double x) {
    return 100 * M_PI * x * std::cos(50 * M_PI * x * x);
  };

  Eigen::MatrixXd results(M_max / dM, 4);
  // Loop over meshes with increasing numbers of cells
  for (int M = dM; M <= M_max; M += dM) {
    // The mesh width
    const double h = 1. / M;

    // Generate the tri-diagonal stiffness matrix for an equidistant
    // mesh on [0, 1] with M cells and p.w. linear Lagrangian finite elements
    // Formulas are explained in Section 2.3 of the lecture document
    Eigen::SparseMatrix<double> A(M + 1, M + 1);
    A.reserve(3 * (M + 1));
    // Fill the diagonal
    for (int i = 0; i < M + 1; ++i) {
      A.coeffRef(i, i) += 2. / h;
    }
    // Fill the off-diagonals
    for (int i = 0; i < M; ++i) {
      A.coeffRef(i, i + 1) += -1. / h;
      A.coeffRef(i + 1, i) += -1. / h;
    }

    // Generate the load vector
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(M + 1);
    for (int i = 0; i < M; ++i) {
      const double a = static_cast<double>(i) / M;
      const Eigen::VectorXd loc_quad_points =
          Eigen::VectorXd::Constant(num_quad_points, a) + h * quad_points;
      const Eigen::VectorXd loc_quad_weights = h * quad_weights;
      // Perform the integration over the cell for both basis functions
      const auto b1 = [&](double x) { return (x - a) / h; };
      const auto b2 = [&](double x) { return 1. - b1(x); };
      for (int k = 0; k < num_quad_points; ++k) {
        rhs[i] += loc_quad_weights[k] * b2(loc_quad_points[k]) *
                  f(loc_quad_points[k]);
        rhs[i + 1] += loc_quad_weights[k] * b1(loc_quad_points[k]) *
                      f(loc_quad_points[k]);
      }
    }

    // Enforce zero dirichlet boundary conditions
    for (long k = 0; k < A.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
        const int row = it.row();
        const int col = it.col();
        if ((row == 0 && col == 0) || (row == M && col == M)) {
          it.valueRef() = 1;
        } else if (row == 0 || row == M || col == 0 || col == M) {
          it.valueRef() = 0;
        }
      }
    }
    // Set the boundary values to zero
    rhs[0] = 0;
    rhs[M] = 0;

    // Solve the resulting linear system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
    const Eigen::VectorXd sol = solver.solve(rhs);

    // Compute the norms and store them in the results matrix
    double norm_max = 0;
    double norm_H1_squared = 0;
    double norm_L2_squared = 0;
    for (int i = 0; i < M; ++i) {
      const double a = static_cast<double>(i) / M;
      const double b = static_cast<double>(i + 1) / M;
      const Eigen::VectorXd loc_quad_points =
          Eigen::VectorXd::Constant(num_quad_points, a) + h * quad_points;
      const Eigen::VectorXd loc_quad_weights = h * quad_weights;

      // The approximate solution on the current cell
      const auto u_h = [&](double x) {
        return sol[i + 1] * (x - a) / h + sol[i] * (1. - (x - a) / h);
      };
      // The gradient of the approximate solution on the current cell
      const auto u_h_grad = [&](double /*x*/) {
        return (sol[i + 1] - sol[i]) / h;
      };
      // The difference of the approximate and the exact solution
      const auto diff = [&](double x) { return u_h(x) - u(x); };
      // The difference in the gradient of the approximate and the exact
      // solution
      const auto diff_grad = [&](double x) { return u_h_grad(x) - u_grad(x); };

      // Compute the max norm by evaluating the functions on a fine grid
      norm_max =
          std::max(norm_max, Eigen::ArrayXd::LinSpaced(10 * M_max / M, a, b)
                                 .unaryExpr(diff)
                                 .abs()
                                 .maxCoeff());
      // Compute the H1 and L2 norms by integrating using a numerical quadrature
      for (int k = 0; k < num_quad_points; ++k) {
        norm_H1_squared += loc_quad_weights[k] * diff_grad(loc_quad_points[k]) *
                           diff_grad(loc_quad_points[k]);
        norm_L2_squared += loc_quad_weights[k] * diff(loc_quad_points[k]) *
                           diff(loc_quad_points[k]);
      }
    }
    results((M / dM) - 1, 0) = M;
    results((M / dM) - 1, 1) = norm_max;
    results((M / dM) - 1, 2) = std::sqrt(norm_H1_squared);
    results((M / dM) - 1, 3) = std::sqrt(norm_L2_squared);
  }

  // Output the resulting errors to a file
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file;
  file.open(output_file);
  file << results.format(CSVFormat);
  file.close();

  return 0;
}
