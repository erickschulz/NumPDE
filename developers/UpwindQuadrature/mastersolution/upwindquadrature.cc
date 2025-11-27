/**
 * @file upwindquadrature.cc
 * @brief NPDE homework template
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "upwindquadrature.h"

namespace UpwindQuadrature {

/**
 * @brief Computes the masses m(p) of all vertices of the mesh
 * @param mesh_p pointer to the mesh.
 * @return Datastructure containing the masses m(p) for all vertices p of the
 * mesh represented by mesh_p.
 */
lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  lf::mesh::utils::CodimMeshDataSet<double> masses(mesh_p, 2, 0.0);
  // compute masses using a cell-based approach.
  for (const lf::mesh::Entity* entity : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_ptr = entity->Geometry();
    double area = lf::geometry::Volume(*geo_ptr);
    for (const lf::mesh::Entity* corner : entity->SubEntities(2)) {
      masses(*corner) += area / 3.0;
    }
  }
  return masses;
}

}  // namespace UpwindQuadrature
