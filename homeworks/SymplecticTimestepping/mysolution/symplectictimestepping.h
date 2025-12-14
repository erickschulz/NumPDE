/**
 * @file symplectictimestepping.h
 * @brief NPDE homework SymplecticTimestepping
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace SymplecticTimestepping {

void sympTimestep(double tau, Eigen::Vector2d& pq_j);

Eigen::Vector2d sympTimesteppingHarmonicOscillatorODE(unsigned int m);

void sympTimesteppingODETest();

Eigen::MatrixXd simulateHamiltonianDynamics(const Eigen::VectorXd& p0,
                                            const Eigen::VectorXd& q0, double T,
                                            unsigned int M);

}  // namespace SymplecticTimestepping
