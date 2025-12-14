/**
 * @file locallaplaceqfe.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "locallaplaceqfe.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>

namespace DebuggingFEM {

Eigen::Matrix<double, 6, 6> LocalLaplaceQFE1::Eval(
    const lf::mesh::Entity& cell) {
  // Query (topological) type of cell/reference element
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Verify that the cell is a triangle
  LF_ASSERT_MSG(ref_el == lf::base::RefEl::kTria(),
                "Implemented for triangles only not for " << ref_el);
  // The final element matrix has size 6x6
  Eigen::Matrix<double, 6, 6> result{};
  // Obtain the vertex coordinates of the triangle
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  // Matrix storing corner coordinates in its columns
  Eigen::Matrix<double, 2, 3> vertices{geo_ptr->Global(ref_el.NodeCoords())};
  // Comopute element matrix for negative Laplacian and lowest-order Lgrangian
  // finite elements as in Remark 2.4.5.9. in the course notes
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  const double area = 0.5 * std::abs(X.determinant());
  // Initialize gradients!
  auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};

  // Returns all gradients of the local shape functions for quadatic Lagrangian
  // finite elements in the columns of a matrix. The gradients are evaluated at
  // a point specified by its reference coordinates.
  auto gradientsLocalShapeFunctions =
      [&grad_bary_coords](
          const Eigen::Vector2d xh) -> Eigen::Matrix<double, 2, 6> {
    Eigen::Matrix<double, 2, 6> gradients;
    // barycentric coordinate functions
    const std::array<double, 3> l{1.0 - xh[0] - xh[1], xh[0], xh[1]};
    gradients.col(0) = grad_bary_coords.col(0) * (4 * l[0] - 1);
    gradients.col(1) = grad_bary_coords.col(1) * (4 * l[1] - 1);
    gradients.col(2) = grad_bary_coords.col(2) * (4 * l[2] - 1);
    gradients.col(3) =
        4 * (grad_bary_coords.col(0) * l[1] + grad_bary_coords.col(1) * l[0]);
    gradients.col(4) =
        4 * (grad_bary_coords.col(1) * l[2] + grad_bary_coords.col(2) * l[1]);
    gradients.col(5) =
        4 * (grad_bary_coords.col(0) * l[2] + grad_bary_coords.col(2) * l[0]);
    return gradients;
  };

  const auto grad_vt_0{gradientsLocalShapeFunctions(Eigen::Vector2d(0, 0))};
  const auto grad_vt_1{gradientsLocalShapeFunctions(Eigen::Vector2d(1, 0))};
  const auto grad_vt_2{gradientsLocalShapeFunctions(Eigen::Vector2d(0, 1))};
  result =
      area / 3.0 *
      (grad_vt_0.transpose() * grad_vt_0 + grad_vt_1.transpose() * grad_vt_1 +
       grad_vt_2.transpose() * grad_vt_2);
  return result;
}

Eigen::Matrix<double, 6, 6> LocalLaplaceQFE2::Eval(
    const lf::mesh::Entity& cell) {
  // Query (topological) type of cell/reference element
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Verify that the cell is a triangle
  LF_ASSERT_MSG(ref_el == lf::base::RefEl::kTria(),
                "Implemented for triangles only not for " << ref_el);
  // The final element matrix has size 6x6
  Eigen::Matrix<double, 6, 6> result{};
  // Obtain the vertex coordinates of the triangle
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  // Matrix storing corner coordinates in its columns
  Eigen::Matrix<double, 2, 3> vertices{geo_ptr->Global(ref_el.NodeCoords())};
  // Comopute element matrix for negative Laplacian and lowest-order Lgrangian
  // finite elements as in Remark 2.4.5.9. in the course notes
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  const double area = 0.5 * std::abs(X.determinant());
  auto grad_bary_coords{X.inverse().block<2, 3>(1, 0)};
  auto L{grad_bary_coords.transpose() * grad_bary_coords};

  // See Example 2.7.5.7 in course notes for derivation of the formulas
  result << 3. * L(0, 0), -L(0, 1), -L(0, 2), 4. * L(0, 1), 0, 4. * L(0, 2),
      -L(0, 1), 3. * L(1, 1), -L(1, 2), 4. * L(0, 1), 4. * L(1, 2), 0, -L(0, 2),
      -L(1, 2), 3. * L(2, 2), 0, 4. * L(2, 1), 4. * L(2, 0), 4. * L(0, 1),
      4. * L(0, 1), 0, 8. * (L(0, 0) + L(0, 1) + L(1, 1)), 8 * L(0, 2),
      8 * L(1, 2), 0, 4. * L(1, 2), 4. * L(2, 1), 8. * L(0, 2),
      8. * (L(1, 1) + L(1, 2) + L(2, 2)), 8 * L(0, 1), 4 * L(0, 2), 0,
      4. * L(2, 0), 8. * L(1, 2), 8. * L(0, 1),
      8. * (L(0, 0) + L(0, 2) + L(2, 2));
  result *= (area / 3.);
  return result;
}

// implementation
Eigen::Matrix<double, 6, 6> LocalLaplaceQFE3::Eval(
    const lf::mesh::Entity& cell) {
  // Obtain the element matrix for piecewise linear Lagrangian FEM by using
  // a built-in class of LehrFEM++
  auto linear_lapl_element_matrix = lf::uscalfe::LinearFELaplaceElementMatrix();
  Eigen::Matrix4d L = linear_lapl_element_matrix.Eval(cell);
  // Variable for returning the final 6x6 element matrix
  Eigen::Matrix<double, 6, 6> result{};

  // The element matrix for quadratic finite elements can be constructed from
  // the element matrix for linear FEM, see Example 2.7.5.7. in the lecture
  // notes.
  result << 3. * L(0, 0), -L(0, 1), -L(0, 2), 4. * L(0, 1), 0, 4. * L(0, 2),
      -L(0, 1), 3. * L(1, 1), -L(1, 2), 4. * L(0, 1), 4. * L(1, 2), 0, -L(0, 2),
      -L(1, 2), 3. * L(2, 2), 0, 4. * L(2, 1), 4. * L(2, 0), 4. * L(0, 1),
      4. * L(0, 1), 0, 8. * (L(0, 0) + L(0, 1) + L(1, 1)), 8 * L(0, 2),
      8 * L(1, 2), 0, 4. * L(1, 2), 4. * L(2, 1), 8. * L(0, 2),
      8. * (L(1, 1) + L(1, 2) + L(2, 2)), 8 * L(0, 1), 4 * L(0, 2), 0,
      4. * L(2, 0), 8. * L(1, 2), 8. * L(0, 1),
      8. * (L(0, 0) + L(0, 2) + L(2, 2));
  // A hideous manipulation introducces an error !
  result(3, 3) *= 1.000001;
  return (result / 3.0);
}

}  // namespace DebuggingFEM
