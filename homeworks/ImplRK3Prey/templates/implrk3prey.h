#ifndef IMPLRK3PREY_H_
#define IMPLRK3PREY_H_

/**
 * @file implrk3prey.h
 * @brief NPDE homework ImplRK3Prey code
 * @author Unknown, Oliver Rietmann
 * @date 29.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>
#include <utility>
#include <vector>

#include "dampnewton.h"

namespace ImplRK3Prey {

// Compute the Kronecker product $C = A \otimes B$
// A is m x n matrix, B is l x k matrix
// return Kronecker product of A and B: dim is m*l x n*k
inline Eigen::MatrixXd kron(const Eigen::MatrixXd& A,
                            const Eigen::MatrixXd& B) {
  Eigen::MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }
  return C;
}

// Implements a Runge-Kutta implicit solver for a given Butcher tableau
// for autonomous ODEs.
/* SAM_LISTING_BEGIN_1 */
class implicitRKIntegrator {
 public:
  // Constructor
  // A is a  matrix containing coefficents of Butcher tableau
  // b is a vector containing coefficients of lower part of Butcher tableau
  implicitRKIntegrator(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
      : A(A), b(b), s(b.size()) {
    assert(A.cols() == A.rows() && "Matrix must be square.");
    assert(A.cols() == b.size() && "Incompatible matrix/vector size.");
  }

  /* Performs the solution of the ODE.
   * Solve an autonomous ODE y' = f(y), y(0) = y0, using an
   * implicit RK scheme given in the Butcher tableau provided in the
   * constructor. Performs N equidistant steps of size h = T / N up to
   * time T with initial condition y0.
   * Returns a vector containing all steps y^n (for each n) inclu. y0 */
  template <class Function, class Jacobian>
  std::vector<Eigen::VectorXd> solve(Function&& f, Jacobian&& Jf, double T,
                                     const Eigen::VectorXd& y0,
                                     unsigned int N) const {
    // Iniz step size
    double h = T / N;

    // Will contain all steps, reserve memory for efficiency
    std::vector<Eigen::VectorXd> res;
    res.reserve(N + 1);

    // Store initial data
    res.push_back(y0);

    // Initialize some memory to store temporary values
    Eigen::VectorXd ytemp1 = y0;
    Eigen::VectorXd ytemp2 = y0;
    // Pointers to swap previous value
    Eigen::VectorXd* yold = &ytemp1;
    Eigen::VectorXd* ynew = &ytemp2;

    // Loop over all fixed steps
    for (unsigned int k = 0; k < N; ++k) {
      // Compute, save and swap next step
      step(f, Jf, h, *yold, *ynew);
      res.push_back(*ynew);
      std::swap(yold, ynew);
    }

    return res;
  }

 private:
  /* Perform a single step of the RK method for the for of the autonomous ODE */
  /* SAM_LISTING_BEGIN_0 */
  template <class Function, class Jacobian>
  void step(Function&& f, Jacobian&& Jf, double h, const Eigen::VectorXd& y0,
            Eigen::VectorXd& y1) const {
    int d = y0.size();
    const Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(d, d);

    // Handle for the function F describing the
    // equation satisfied by the stages g
    auto F = [&y0, h, d, this, &f, &eye](Eigen::VectorXd gv) {
      Eigen::VectorXd Fv = gv;
      for (int j = 0; j < s; j++) {
        Fv = Fv - h * kron(A.col(j), eye) * f(y0 + gv.segment(j * d, d));
      }
      return Fv;
    };

    // Handle for the Jacobian of F.
    auto JF = [&y0, h, d, &Jf, this, &eye](Eigen::VectorXd gv) {
      Eigen::MatrixXd DF(s * d, s * d);
      for (int j = 0; j < s; j++) {
        DF.block(0, j * d, s * d, d) =
            kron(A.col(j), eye) * Jf(y0 + gv.segment(j * d, d));
      }
      DF = Eigen::MatrixXd::Identity(s * d, s * d) - h * DF;
      return DF;
    };

    //====================
    // Your code goes here
    // Use the function dampnewton(...) in dampnewton.h
    // to solve the non-linear system of equations.
    //====================
  }
  /* SAM_LISTING_END_0 */
  //<! Matrix A in Butcher scheme
  const Eigen::MatrixXd A;
  //<! Vector b in Butcher scheme
  const Eigen::VectorXd b;
  //<! Size of Butcher matrix and vector A and b
  unsigned int s;
};
/* SAM_LISTING_END_1 */

}  // namespace ImplRK3Prey

#endif  // #define IMPLRK3PREY_H_
