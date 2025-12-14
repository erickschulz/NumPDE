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
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

}  // namespace DebuggingFEM

#endif
