#ifndef IMPLICITRKINTEGRATOR_H_
#define IMPLICITRKINTEGRATOR_H_

/**
 * @file implicitrkintegrator.h
 * @brief NPDE homework CrossProd code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>
#include <utility>
#include <vector>

#include "dampnewton.h"

namespace CrossProd {

// Compute the Kronecker product A and B.
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

// Implements a Runge-Kutta implicit solver for a
// given Butcher tableau for autonomous ODEs.
class implicitRKIntegrator {
 public:
  // Constructor for the implicit RK method.
  implicitRKIntegrator(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
      : A(A), b(b), s(b.size()) {
    assert(A.cols() == A.rows() && "Matrix must be square.");
    assert(A.cols() == b.size() && "Incompatible matrix/vector size.");
  }

  /* Perform the solution of the ODE.
  Solve an autonomous ODE y' = f(y), y(0) = y0, using an
  implicit RK scheme given in the Butcher tableau provided in the
  constructor. Performs N equidistant steps upto time T
  with initial data y0. */
  template <class Function, class Jacobian>
  std::vector<Eigen::VectorXd> solve(Function&& f, Jacobian&& Jf, double T,
                                     const Eigen::VectorXd& y0,
                                     unsigned int M) const {
    // Iniz step size
    double h = T / M;

    // Will contain all steps, reserve memory for efficiency
    std::vector<Eigen::VectorXd> res;
    res.reserve(M + 1);

    // Store initial data
    res.push_back(y0);

    // Initialize some memory to store temporary values
    Eigen::VectorXd ytemp1 = y0;
    Eigen::VectorXd ytemp2 = y0;
    // Pointers to swap previous value
    Eigen::VectorXd* yold = &ytemp1;
    Eigen::VectorXd* ynew = &ytemp2;

    // Loop over all fixed steps
    for (unsigned int k = 0; k < M; ++k) {
      // Compute, save and swap next step
      step(std::forward<Function>(f), std::forward<Jacobian>(Jf), h, *yold,
           *ynew);
      res.push_back(*ynew);
      std::swap(yold, ynew);
    }

    return res;
  }

 private:
  // Perform a single step of the RK method for the sol. of the autonomous ODE
  // Compute a single explicit RK step y^{n+1} = y_n + \sum ...
  // starting from value y0 and storing next value in y1
  /* SAM_LISTING_BEGIN_0 */
  template <class Function, class Jacobian>
  void step(Function&& f, Jacobian&& Jf, double h, const Eigen::VectorXd& y0,
            Eigen::VectorXd& y1) const {
    int d = y0.size();
    const Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(d, d);

    // Handle for the function F describing the
    // equation satisfied by the stages g
    auto F = [&y0, h, d, this, &f,
              &eye](Eigen::VectorXd gv) -> Eigen::VectorXd {
      Eigen::VectorXd Fv = gv;
      for (int j = 0; j < s; j++) {
        Fv = Fv - h * kron(A.col(j), eye) * f(y0 + gv.segment(j * d, d));
      }
      return Fv;
    };

    // Handle for the Jacobian of F.
    auto JF = [&y0, h, d, &Jf, this,
               &eye](Eigen::VectorXd gv) -> Eigen::MatrixXd {
      Eigen::MatrixXd DF(s * d, s * d);
      for (int j = 0; j < s; j++) {
        DF.block(0, j * d, s * d, d) =
            kron(A.col(j), eye) * Jf(y0 + gv.segment(j * d, d));
      }
      DF = Eigen::MatrixXd::Identity(s * d, s * d) - h * DF;
      return DF;
    };

    // Obtain stages with damped Newton method
    Eigen::VectorXd gv = Eigen::VectorXd::Zero(s * d);
    dampnewton(F, JF, gv);

    // Calculate y1
    Eigen::MatrixXd K(d, s);
    for (int j = 0; j < s; j++) K.col(j) = f(y0 + gv.segment(j * d, d));
    y1 = y0 + h * K * b;
  }
  /* SAM_LISTING_END_0 */
  //<! Matrix A in Butcher scheme
  const Eigen::MatrixXd A;
  //<! Vector b in Butcher scheme
  const Eigen::VectorXd b;
  //<! Size of Butcher matrix and vector A and b
  unsigned int s;
};

}  // namespace CrossProd

#endif  // #define IMPLICITRKINTEGRATOR_H_
