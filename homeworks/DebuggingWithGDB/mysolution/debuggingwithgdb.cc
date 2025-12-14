/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "debuggingwithgdb.h"

#include <cstdlib>

#include "systemcall.h"

namespace DebuggingWithGDB {
/* SAM_LISTING_BEGIN_1 */
void ReadAndOutputMesh(const char* filename) {
  if (filename != nullptr) {
    // Build full path to the mesh file
    auto gmshfile_path = filename;
    // Load the mesh from a file produced by Gmsh
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), gmshfile_path);
    // Obtain pointer to read mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
    const lf::mesh::Mesh& mesh{*mesh_p};
    // Run through all entities of the mesh and print entity information
    for (int codim = 0; codim <= mesh.DimMesh(); ++codim) {
      int cnt = 0;
      Eigen::VectorXd c{Eigen::VectorXd::Zero(mesh.DimWorld())};
      for (const lf::mesh::Entity* entity : mesh.Entities(codim)) {
        // Number of vertices
        const int no_vertices = (entity->RefEl()).NumNodes();
        // Obtain "convex hull" of an entity
        const Eigen::MatrixXd vertices{
            lf::geometry::Corners(*entity->Geometry())};
        // Compute center of gravity
        c += vertices.rowwise().sum() / no_vertices;
        cnt++;
      }
      if (cnt != mesh.NumEntities(codim)) {
        std::cerr << "Count mismatch for entities of codim = " << codim << ": "
                  << cnt << " <-> " << mesh.NumEntities(codim) << std::endl;
      }
      c /= cnt;
      std::cout << "Center of codim-" << codim
                << " entities = " << c.transpose() << std::endl;
    }
    // Wite mesh data to file for visualization with Python script
    lf::io::writeMatplotlib(mesh, "ljoint.csv");
    std::cout << "Wrote ljoint.csv" << std::endl;
    systemcall::execute("python3 scripts/plot_mesh.py ljoint.csv mesh.eps");
  }
}
/* SAM_LISTING_END_1 */

}  // namespace DebuggingWithGDB
