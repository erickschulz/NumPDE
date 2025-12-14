#include "../simplelinearfiniteelements.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <tuple>

#include "../tria_mesh_2D.h"

namespace SimpleLinearFiniteElements::test {

constexpr double pi = 3.1415926535897;
constexpr char meshfile1[] = "meshes/Square1.txt";
constexpr char meshfile3[] = "meshes/Square3.txt";

/**
 * @brief test ElmentMatrix_Mass_LFE implementation
 */
TEST(SimpleLinearFiniteElements, ElementMatrix_Mass_LFE) {
  // check the produced matrix for a fairly standard triangle
  Eigen::Matrix<double, 2, 3> test;
  test << 0, 1, 0, 0, 0, 1;
  Eigen::Matrix3d M;
  M = ElementMatrix_Mass_LFE(test);

  Eigen::Matrix3d ref_M;
  ref_M << 0.0833333, 0.0416667, 0.0416667, 0.0416667, 0.0833333, 0.0416667,
      0.0416667, 0.0416667, 0.0833333;

  double tol = 1e-8;
  ASSERT_NEAR(ref_M.norm(), M.norm(), tol);
}

/**
 * @brief test L2Error implementation
 */
TEST(SimpleLinearFiniteElements, L2Error) {
  // read coarsest mesh
  TriaMesh2D square_mesh(meshfile1);
  // exact solution
  auto uExact = [](const Eigen::Vector2d& x) {
    return std::cos(2 * pi * x(0)) * std::cos(2 * pi * x(1));
  };
  // source function
  auto f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // assemble galerkin matrix and load vector
  Eigen::SparseMatrix<double> A =
      assembleGalMatLFE(square_mesh, getElementMatrix);
  Eigen::VectorXd L = assemLoad_LFE(square_mesh, localLoadLFE, f);

  // solve linear system of equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Eigen::VectorXd U = solver.solve(L);

  // compare to expected error
  double error = L2Error(square_mesh, U, uExact);
  ASSERT_NEAR(error, 0.232547, 0.1);
}

/**
 * @brief test H1Serror implementation
 */
TEST(SimpleLinearFiniteElements, H1Serror) {
  // read coarsest mesh
  TriaMesh2D square_mesh(meshfile3);
  // exact gradient
  auto gradUExact = [](const Eigen::Vector2d& x) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x(0)) * std::cos(2 * pi * x(1)),
        -2 * pi * std::cos(2 * pi * x(0)) * std::sin(2 * pi * x(1));
    return gradient;
  };
  // source function
  auto f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };
  // compute galerkin matrix and load vector
  Eigen::SparseMatrix<double> A =
      assembleGalMatLFE(square_mesh, getElementMatrix);
  Eigen::VectorXd L = assemLoad_LFE(square_mesh, localLoadLFE, f);

  // solve linear system of equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Eigen::VectorXd U = solver.solve(L);

  // compare to expected error
  // high tolerance as the procedure is affected by many rounding errors
  double error = H1Serror(square_mesh, U, gradUExact);
  ASSERT_NEAR(error, 1.32457, 0.1);
}

/**
 * @brief test solve implementation
 */
TEST(SimpleLinearFiniteElements, Solve) {
  TriaMesh2D square_mesh(meshfile1);
  std::tuple<Eigen::VectorXd, double, double> solution = Solve(square_mesh);

  Eigen::VectorXd solution_ref(25);
  solution_ref << 0.625091201454706, 1.82631947928602, 0.625091201454444,
      1.82631947928602, -1.25162595870465, -1.25162595870465, 1.22710538961478,
      -1.25162595870452, -1.25162595870452, -0.205607688673624,
      -0.205607688673623, 0.0791862852828167, 0.180931017968087,
      0.18093101796835, -0.205607688673888, -0.205607688673888,
      0.180931017968087, 0.079186285282817, 0.180931017968351,
      -0.104406878994863, -0.0124922491782952, -0.012492249178295,
      -0.0124922491780322, -0.0124922491780322, -0.104406878994864;
  double L2Error_ref = 0.232547;
  double H1Serror_ref = 4.95432;

  ASSERT_EQ(std::get<0>(solution).size(), solution_ref.size());
  double tol = 1.0e-5;
  ASSERT_NEAR(0.0,
              (std::get<0>(solution) - solution_ref).lpNorm<Eigen::Infinity>(),
              tol);
  ASSERT_NEAR(0.0, std::get<1>(solution) - L2Error_ref, tol);
  ASSERT_NEAR(0.0, std::get<2>(solution) - H1Serror_ref, tol);
}

}  // namespace SimpleLinearFiniteElements::test
