/**
 * @file neumanndatarecovery.h
 * @brief NPDE homework NeumannDataRecovery code
 * @author R. Hiptmair
 * @date July 14, 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef NDR_H_
#define NDR_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace NeumannDataRecovery {
/** @brief Computation of gradients of barycentric coordinate functions
 *
 * @param vertices 2x3 matrix whose columns contain the vertex coordinates of
 * the triangle
 * @return gradients of barycentric coordinate functions stored in the columns
 * of a 2x3 matrix.
 *
 * For explanations see Code 2.4.5.11 in the NumPDE lecture document.
 */
Eigen::Matrix<double, 2, 3> GradsBaryCoords(
    Eigen::Matrix<double, 2, 3> vertices);

/** @brief Computation of exterior unit normals for flat triangle/quadrilateral
 *
 * @param geo reference to the goemetry of the entity of co-dimension 0 for
 * which the exterior unit normals are requested.
 * @return 2xn-matrix, n = no. of vertices, whose columns contain the unit
 * normals.
 *
 * @note the function always returns a 2x4 matrix. In the case of a triangle the
 * last column is just not initialized.
 *
 * Algorithm: For a triangle the gradients of the barycentric coordinate
 * functions provided the directions of the (negative) exterior normals. A
 * quadrilateral is decomposed into two triangles, whose exterior unit normals
 * are computed subsequently.
 */
Eigen::Matrix<double, 2, 4> exteriorUnitNormals(
    const lf::geometry::Geometry& geo);

/** @brief Compute the piecewise constant Neumann data directly from a piecewise
 * linear finite-element solution
 *
 * @param fe_space reference to lowest order Lagrangian finite element space on
 * a triangular mesh
 * @param mu basis expansion coefficient vector (covering all basis functions)
 * @return constant value of Neumann trace for every edge on the boundary, zero
 * for interior edges.
 *
 */
lf::mesh::utils::CodimMeshDataSet<double> getNeumannData(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& mu);

}  // namespace NeumannDataRecovery

#endif
