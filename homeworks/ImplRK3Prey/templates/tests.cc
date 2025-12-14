#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    // Dimension of state space
    d = 2;

    // Variables for testing (iterations, final times)
    N = {128, 256, 512, 1024, 2048};
    T = {2.0, 5.0, 10.0, 20.0, 50.0};

    // Initial value for model
    y0.resize(d);
    y0 << 100, 5;

    // Definition of coefficients in Butcher scheme
    s = 2;

    A.resize(s, s);
    b.resize(s);
    A << 5. / 12., -1. / 12., 3. / 4., 1. / 4.;
    b << 3. / 4., 1. / 4.;

    // Coefficients for predator/prey model
    alpha1 = 3.;
    alpha2 = 2.;
    beta1 = 0.1;
    beta2 = 0.1;

    f = [this](const VectorXd& y) {
      auto temp = y;
      temp(0) *= alpha1 - beta1 * y(1);
      temp(1) *= -alpha2 + beta2 * y(0);
      return temp;
    };

    Jf = [this](const VectorXd& y) {
      MatrixXd temp(2, 2);
      temp << alpha1 - beta1 * y(1), -beta1 * y(0), beta2 * y(1),
          -alpha2 + beta2 * y(0);
      return temp;
    };
  }

  unsigned int s;
  unsigned int d;

  std::vector<unsigned int> N;
  std::vector<double> T;

  VectorXd y0;
  MatrixXd A;
  VectorXd b;

  double alpha1;
  double alpha2;
  double beta1;
  double beta2;

  std::function<MatrixXd(const VectorXd&)> Jf;
  std::function<VectorXd(const VectorXd&)> f;
};

TestData data;

TEST_SUITE("ImplRK Prey") {
  TEST_CASE("std::vector<VectorXd> solve" *
            doctest::description("RK solving")) {
    // Instantiate objects
    implicit_RKIntegrator RK_sol(data.A, data.b);
    implicit_RKIntegrator_TEST RK_stud(data.A, data.b);

    // Test with various values
    for (int i = 0; i < data.T.size(); i++) {
      for (int j = 0; j < data.N.size(); j++) {
        std::vector<VectorXd> sol_vec =
            RK_sol.solve(data.f, data.Jf, data.T[i], data.y0, data.N[j]);
        std::vector<VectorXd> stud_vec =
            RK_stud.solve_TEST(data.f, data.Jf, data.T[i], data.y0, data.N[j]);

        VectorXd sol = sol_vec.back();
        VectorXd stud = stud_vec.back();

        CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
      }
    }

    MESSAGE("Run the program to see a convergence study.");
  }

  TEST_CASE("MatrixXd kron" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("implicit_RKIntegrator" * doctest::description("skipped") *
            doctest::skip()) {}
}
