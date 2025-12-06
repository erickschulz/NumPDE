/**
 * @file nonlinschroedingerequation.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "nonlinschroedingerequation.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <cmath>

namespace NonLinSchroedingerEquation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix3d MassElementMatrixProvider::Eval(const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  Eigen::Matrix3d element_matrix;
  //====================
  // Your code goes here
  //====================
  return element_matrix;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double Norm(const Eigen::VectorXcd& mu, const Eigen::SparseMatrix<double>& D) {
  //====================
  // Your code goes here
  // Replace this dummy value by the approximate
  // norm of the mesh function assiciated with mu
  return 1.0;
  //====================
}

double KineticEnergy(const Eigen::VectorXcd& mu,
                     const Eigen::SparseMatrix<double>& A) {
  //====================
  // Your code goes here
  // Replace this dummy value by the kinetic energy
  // of the mesh function assiciated with mu
  return 0.0;
  //====================
}

double InteractionEnergy(const Eigen::VectorXcd& mu,
                         const Eigen::SparseMatrix<double>& D) {
  //====================
  // Your code goes here
  // Replace this dummy value by the interaction
  // energy of the mesh function assiciated with mu
  return 0.0;
  //====================
}
/* SAM_LISTING_END_2 */

}  // namespace NonLinSchroedingerEquation
