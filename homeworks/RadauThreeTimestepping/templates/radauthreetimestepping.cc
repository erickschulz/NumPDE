/**
 * @file radauthreetimestepping.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz, edited by Oliver Rietmann
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimestepping.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>

namespace RadauThreeTimestepping {

/**
 * @brief Implementation of the right hand side (time dependent) source vector
 * for the parabolic heat equation
 * @param dofh A reference to the DOFHandler
 * @param time The time at which to evaluate the source vector
 * @returns The source vector at time `time`
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler& dofh,
                                    double time) {
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Right-hand side vector has to be set to zero initially
  Eigen::VectorXd phi(N_dofs);
  //====================
  // Your code goes here
  //====================
  return phi;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Heat evolution solver: the solver obtains the
 * discrete evolution operator from the Radau3MOLTimestepper class and
 * repeatedly iterates its applicaiton starting from the initial condition
 * @param dofh The DOFHandler object
 * @param m is total number of steps until final time final_time (double)
 * @param final_time The duration for which to solve the PDE
 * @returns The solution at the final timestep
 */
/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler& dofh,
                                   unsigned int m, double final_time) {
  Eigen::VectorXd discrete_heat_sol(dofh.NumDofs());
  //====================
  // Your code goes here
  //====================
  return discrete_heat_sol;
}
/* SAM_LISTING_END_6 */

/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity& tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  //====================
  // Your code goes here
  //====================
  return elMat;  // return the local mass element matrix
}

/* Implementing constructor of class Radau3MOLTimestepper */
/* SAM_LISTING_BEGIN_4 */
Radau3MOLTimestepper::Radau3MOLTimestepper(const lf::assemble::DofHandler& dofh)
    : dofh_(dofh) {
  //====================
  // Your code goes here
  // Add any additional members you need in the header file
  //====================
}
/* SAM_LISTING_END_4 */

/* Implementation of Radau3MOLTimestepper member functions */
// The function discreteEvolutionOperator() returns the discretized evolution
// operator as obtained from the Runge-Kutta Radau IIA 2-stages method using the
// Butcher table as stored in the Radau3MOLTimestepper class
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd Radau3MOLTimestepper::discreteEvolutionOperator(
    double time, double tau, const Eigen::VectorXd& mu) const {
  Eigen::VectorXd discrete_evolution_operator(dofh_.NumDofs());
  //====================
  // Your code goes here
  //====================
  return discrete_evolution_operator;
}
/* SAM_LISTING_END_5 */

}  // namespace RadauThreeTimestepping
