/**
 * @file potentialflow.h
 * @brief NPDE homework Potential Flow code
 * @author R. Hiptmair
 * @date January 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef POTFLOW_H_
#define POTFLOW_H_
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace PotentialFlow {
// Initialization of Galerkin matrix
extern Eigen::SparseMatrix<double> initializeA(unsigned int M);

// Initialization of r.h.s. vector based on Dirichlet data passed in g
extern Eigen::VectorXd initializeRHSVector(
    const std::function<double(double, double)>& g, unsigned int M);

}  // namespace PotentialFlow

#endif
