/**
 * @file neumanndatarecovery.cc
 * @brief NPDE homework NeumannDataRecovery code
 * @author R. Hiptmair
 * @date July 14, 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#define _USE_MATH_DEFINES
#include "neumanndatarecovery.h"

#include <math.h>

#include <cmath>

namespace NeumannDataRecovery {
Eigen::Matrix<double, 2, 3> GradsBaryCoords(
    Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9 in the lecture document
  Eigen::Matrix<double, 3, 3> X;
  // See (2.4.5.10) in lecture document: first column of X contains all 1s, the
  // other two columns the first and second coordinates of the vertex coordinate
  // vectors
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // as columns of a 2x3 matrix containing the \beta-coefficients in (2.4.5.10).
  return X.inverse().block<2, 3>(1, 0);
}

Eigen::Matrix<double, 2, 4> exteriorUnitNormals(
    const lf::geometry::Geometry& geo) {
  // Only available for flat triangles and quadrilaterals
  LF_ASSERT_MSG(geo.DimGlobal() == 2, "Only implemented for 2D meshes");
  // Return variable
  Eigen::Matrix<double, 2, 4> unit_normals;
  const lf::base::RefEl ref_el{geo.RefEl()};
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      const Eigen::Matrix<double, 2, 3> vertex_coords{
          lf::geometry::Corners(geo)};
      const Eigen::Matrix<double, 2, 3> grad_bc{GradsBaryCoords(vertex_coords)};
      // Number unit normals according to the numbering of edges
      unit_normals.col(0) = -grad_bc.col(2).normalized();
      unit_normals.col(1) = -grad_bc.col(0).normalized();
      unit_normals.col(2) = -grad_bc.col(1).normalized();
      unit_normals.col(3) = Eigen::Vector2d::Zero();
      break;
    }
    case lf::base::RefEl::kQuad(): {
      const Eigen::Matrix<double, 2, 4> vc{lf::geometry::Corners(geo)};
      // Split into two triangles
      const Eigen::Matrix<double, 2, 3> T1 = vc.block<2, 3>(0, 0);
      const Eigen::Matrix<double, 2, 3> T2 =
          (Eigen::Matrix<double, 2, 3>() << vc.col(0), vc.block<2, 2>(0, 2))
              .finished();
      const Eigen::Matrix<double, 2, 3> gbc1{GradsBaryCoords(T1)};
      const Eigen::Matrix<double, 2, 3> gbc2{GradsBaryCoords(T2)};
      // Number unit normals according to the numbering of edges
      unit_normals.col(0) = -gbc1.col(2).normalized();
      unit_normals.col(1) = -gbc1.col(0).normalized();
      unit_normals.col(2) = -gbc2.col(0).normalized();
      unit_normals.col(3) = -gbc2.col(1).normalized();
      break;
    }
    default: {
      LF_VERIFY_MSG(false, "Unsupported cell type " << ref_el);
      break;
    }
  }
  return unit_normals;
}

/* SAM_LISTING_BEGIN_1 */
lf::mesh::utils::CodimMeshDataSet<double> getNeumannData(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& mu) {
  lf::mesh::utils::CodimMeshDataSet<double> edge_vals(fe_space->Mesh(), 1, 0.0);
  // ========================================
  // Your code here
  // ========================================
  return edge_vals;
}
/* SAM_LISTING_END_1 */

}  // namespace NeumannDataRecovery
