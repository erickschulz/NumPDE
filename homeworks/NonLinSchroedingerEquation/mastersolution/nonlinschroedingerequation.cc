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
  double area = lf::geometry::Volume(*(cell.Geometry()));
  element_matrix = area / 3.0 * Eigen::Matrix3d::Identity();
  return element_matrix;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double Norm(const Eigen::VectorXcd& mu, const Eigen::SparseMatrix<double>& D) {
  return std::sqrt(mu.dot(D * mu).real());
}

double KineticEnergy(const Eigen::VectorXcd& mu,
                     const Eigen::SparseMatrix<double>& A) {
  return 0.5 * mu.dot(A * mu).real();
}

double InteractionEnergy(const Eigen::VectorXcd& mu,
                         const Eigen::SparseMatrix<double>& D) {
  Eigen::VectorXd mu_abs2 = mu.cwiseAbs2();
  return 0.25 * mu_abs2.dot(D * mu_abs2);
}
/* SAM_LISTING_END_2 */

}  // namespace NonLinSchroedingerEquation
