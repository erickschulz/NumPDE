#ifndef ORDNOTALL_H_
#define ORDNOTALL_H_

#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <vector>

#include "rkintegrator.h"

namespace OrdNotAll {

/*!
 * \brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 * \tparam Function Type for r.h.s function f.
 * \param f The r.h.s function for the ODE.
 * \param T Final time.
 * \param y0 Initial data.
 * \param A Butcher matrix $A$.
 * \param b Butcher vector $b$.
 */
template <class Function>
void testCvgRKSSM(const Function& f, double T, double y0,
                  const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

void cmpCvgRKSSM();

}  // namespace OrdNotAll

#endif
