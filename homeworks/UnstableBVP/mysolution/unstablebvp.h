/**
 * @file
 * @brief Solution of source-free heat equation and computation of H1
 *  	  seminorms on different triangular meshes and refinement levels
 * @author Julien Gacon, Amélie Loher
 * @date   March 2019
 * @copyright MIT License
 */

#include <memory>
#include <string>
// Lehrfempp
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

namespace UnstableBVP {

/** @brief Create a mesh hierarchy with `reflevels` levels for a mesh
 *        in shape of a triangle. The triangle can be above x2 = 0 ("top"),
 *        below ("bottom") or intersecting it (any other type)
 * @param reflevels Number of refinement levels for mesh
 * @param mesh_type Where triangle is located, see description
 * @return A shared pointer to a lf::refinement::MeshHierarchy object
 */
std::shared_ptr<lf::refinement::MeshHierarchy> createMeshHierarchy(
    const int reflevels, const std::string& mesh_type = "top");

/** @brief Solve source-free diffusion PDE with a potentially non-continuous
 *        boundary condition, depending on the x2 variable, and compute H1
 *        seminorm of the solution
 * @param mesh_p The mesh, a shared pointer to a `lf::mesh::Mesh object`
 * We have to pass a pointer, because it is required by the constructor of
 * a `lf::uscalfe::FeSpaceLagrangeO1` object.
 * @return H1 seminorm of the solution
 */
double solveTemperatureDistribution(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

}  // namespace UnstableBVP
