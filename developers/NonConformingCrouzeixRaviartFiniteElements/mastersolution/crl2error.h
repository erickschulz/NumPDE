/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   18.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_COMPUTE_CR_L2_ERROR_H
#define NUMPDE_COMPUTE_CR_L2_ERROR_H

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include "crfespace.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTION>
double computeCRL2Error(std::shared_ptr<CRFeSpace> fe_space,
                        const Eigen::VectorXd& mu, FUNCTION&& u) {
  double l2_error = 0.;

// TODO: task 2-14.w)
#if SOLUTION
  // Obtain local-to-global map and current mesh object
  const lf::assemble::DofHandler& dof_handler{fe_space->LocGlobMap()};
  auto mesh_ptr = fe_space->Mesh();

  // Loop over all cells of the mesh (entities of co-dimension 0)
  for (const lf::mesh::Entity* cell : mesh_ptr->Entities(0)) {
    const lf::assemble::size_type num_nodes = cell->RefEl().NumNodes();
    LF_ASSERT_MSG(num_nodes == 3, "Only meaningful for triangles!");
    // Obtain pointer to shape information for cell
    const lf::geometry::Geometry* cell_geom{cell->Geometry()};
    // 2x3- matrix with corner coordinates in its columns
    const Eigen::MatrixXd vertices{lf::geometry::Corners(*cell_geom)};
    // clang-format off
    // 2x3-matrix of midpoint coordinates
    Eigen::MatrixXd midpoints{vertices *
      (Eigen::Matrix<double, 3, 3>(3,3) <<
				 0.5, 0.0, 0.5,
				 0.5, 0.5, 0.0,
				 0.0, 0.5, 0.5)
		     .finished()};
    // clang-format on
    // Obtain global indices of local shape functions
    std::span<const lf::assemble::gdof_idx_t> cell_dof_idx(
        dof_handler.GlobalDofIndices(*cell));

    // Sum contributions of quadrature nodes
    double local_sum = 0.;
    // The CR interpolation nodes are the midpoints and so the exact
    // solution needs to be evaluated at the same points
    for (int loc_idx = 0; loc_idx < num_nodes; ++loc_idx) {
      local_sum +=
          std::pow(mu[cell_dof_idx[loc_idx]] - u(midpoints.col(loc_idx)), 2);
    }
    // Sum cell contributions
    l2_error += lf::geometry::Volume(*cell_geom) * (local_sum / num_nodes);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return std::sqrt(l2_error);
}
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_COMPUTE_CR_L2_ERROR_H
