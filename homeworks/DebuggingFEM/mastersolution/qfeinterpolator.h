/**
 * @file qfeinterpolator.h
 * @brief NPDE homework DebuggingFEM code
 * @author Simon Meierhans & Ralf Hiptmair
 * @date 27/03/2019 & 27.03.2024
 * @copyright Developed at ETH Zurich
 */

#ifndef NPDECODES_DEBUGGINGFEM_QFEINTERPOLATOR_H_
#define NPDECODES_DEBUGGINGFEM_QFEINTERPOLATOR_H_

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace DebuggingFEM {

using size_type = lf::base::size_type;

/**
 * @brief get the global coordinates of a interpolation point in a cell
 * @param idx local index
 *        cell the cell to be used
 * @returns The global coordinate of the i-th interpolation node in the given
 * cell
 */
Eigen::Vector2d globalCoordinate(int idx, const lf::mesh::Entity& cell);

/**
 * @brief interpolate function over a second order lagrangian finite element
 * space
 * @param dofh dof-handler
 *        f function to interpolate
 * @returns A vector containing thebasis function coefficients
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
Eigen::VectorXd interpolateOntoQuadFE(const lf::assemble::DofHandler& dofh,
                                      FUNCTOR&& f) {
  // Obtain a pointer to the mesh object
  auto mesh = dofh.Mesh();
  const size_type N_dofs(dofh.NumDofs());
  // variable for returning result vector
  Eigen::VectorXd result = Eigen::VectorXd::Zero(N_dofs);
  // Reference coordinates of the interpolation nodes of the triangle
  Eigen::Matrix<double, 2, 6> refnodes(2, 6);
  refnodes << 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5;
  // Loop over the cells of the mesh (codim-0 entities)
  for (const lf::mesh::Entity* cell : mesh->Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Implemented for triangles only");
    // Fetch pointer to asscoiated geometry object
    const lf::geometry::Geometry* geom = cell->Geometry();
    // Obtain actual coordinates of interpolation nodes
    Eigen::MatrixXd nodes{geom->Global(refnodes)};
    // get local to global index map for the current cell
    auto glob_ind = dofh.GlobalDofIndices(*cell);
    // Loop over local interpolation nodes
    for (int i = 0; i < 6; i++) {
      // update the result vector
      result(glob_ind[i]) = f(nodes.col(i));
    }
  }
  return result;
}
/* SAM_LISTING_END_1 */

}  // namespace DebuggingFEM

#endif
