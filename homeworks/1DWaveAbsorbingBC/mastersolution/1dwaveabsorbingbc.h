/**
 * @file 1dwaveabsorbingbc.h
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <utility>

namespace WaveAbsorbingBC1D {

/**
 * @brief Computes the solution u(x,t) of the wave equation with paramter c on
 * [0,1]x[0,T] using the leapfrog scheme.
 * @param c speed of propagation
 * @param T endtime
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param m number of timesteps
 * @return Matrix R of shape (m+1) times (N+1) s.t. R(i,j) = u(x_j, t_i)
 */
Eigen::MatrixXd waveLeapfrogABC(double c, double T, unsigned int N,
                                unsigned int m);

/**
 * @brief Computes the potential and kinetic energy of the above solution u
 * @param full_solution matrix R from above
 * @param c speed of propagation
 * @param tau timestep-size (in the notation above: tau = T / m)
 * @return Pair of potential and kinetic energy (in this order) at times
 * 0=t_0<...<t_{m}=T
 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> computeEnergies(
    const Eigen::MatrixXd& full_solution, double c, double tau);

}  // namespace WaveAbsorbingBC1D
