/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;
#if SOLUTION
  // iterate over all cells (co-dimension = 0)
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    volume += lf::geometry::Volume(*geo_p);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
#if SOLUTION
  // Obtain an array of boolean flags for the vertices of the mesh: 'true'
  // indicates that the vertex lies on the boundary.
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // iterate over all edges (co-dimension = 1)
  for (const lf::mesh::Entity* cell : mesh_p->Entities(1)) {
    // check if edge is part of the boundary
    if (bd_flags(*cell)) {
      const lf::geometry::Geometry* geo_p = cell->Geometry();
      length += lf::geometry::Volume(*geo_p);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string filename) {
  double volume, length;

#if SOLUTION
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), filename);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  // call the functions we already implemented
  length = LengthOfBoundary::lengthOfBoundary(mesh_p);
  volume = LengthOfBoundary::volumeOfDomain(mesh_p);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
