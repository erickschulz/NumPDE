/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearloadvector.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <functional>

namespace ElementMatrixComputation {

namespace {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector4d computeLoadVector(
    const Eigen::MatrixXd& vertices,
    std::function<double(const Eigen::Vector2d&)> f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();
  // Vector for returning element vector
  Eigen::Vector4d elem_vec = Eigen::Vector4d::Zero();

  //====================
  // Your code goes here
  //====================

  return elem_vec;
}
/* SAM_LISTING_END_1 */

}  // namespace

Eigen::Vector4d MyLinearLoadVector::Eval(const lf::mesh::Entity& cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

}  // namespace ElementMatrixComputation
