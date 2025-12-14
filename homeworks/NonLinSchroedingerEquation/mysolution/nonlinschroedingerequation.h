#ifndef NONLINSCHROEDINGEREQUATION_H_
#define NONLINSCHROEDINGEREQUATION_H_

/**
 * @file nonlinschroedingerequation.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace NonLinSchroedingerEquation {

/** @brief Class providing the local element matrix needed for assembly.
 *  It satisfies the LehrFEM++ concept EntityMatrixProvider.
 */
class MassElementMatrixProvider {
 public:
  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity& cell) { return true; }
  /** @brief routine for the computation of element matrices
   *  @param cell reference to the triangular cell for
   *         which the element matrix should be computed.
   *  @return element matrix
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity& cell);
};

/** @brief Computes the $L^2(\mathbb{R^2;\mathbb{C}})$-norm, approximated by 2D
 * trapezoidal rule.
 *  @param mu vector of length $N$ containing nodal values
 *  @param D real mass matrix of shape $N \times N$
 *  @return $L^2$-norm of the complex-valued mesh function represented by mu
 */
double Norm(const Eigen::VectorXcd& mu, const Eigen::SparseMatrix<double>& D);

/** @brief Computes the kinetic energy.
 *  @param mu vector of length $N$ containing nodal values
 *  @param A Galerkin matrix of $-\Delta$ (stiffness matrix) of shape $N \times
 * N$
 *  @return kinetic energy of the complex-valued mesh function associated with
 * mu
 */
double KineticEnergy(const Eigen::VectorXcd& mu,
                     const Eigen::SparseMatrix<double>& A);

/** @brief Computes the interaction energy (i.e. energy associated with the
 * non-linear term).
 *  @param mu vector of length $N$ containing nodal values
 *  @param D real mass matrix of shape $N \times N$
 *  @return interaction energy of the complex-valued mesh function associated
 * with mu
 */
double InteractionEnergy(const Eigen::VectorXcd& mu,
                         const Eigen::SparseMatrix<double>& D);

}  // namespace NonLinSchroedingerEquation

#endif  // NONLINSCHROEDINGEREQUATION_H_
