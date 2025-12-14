/**
 * @file gausslobattoparabolic_main.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "gausslobattoparabolic.h"

constexpr double PI = 3.14159265358979323846;

int main() {
  // Load the mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/hex_hybrid.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Solve the PDE
  double T = 1.0;
  unsigned int M = 100;
  auto g = [](double t) { return t < 1.0 ? std::sin(0.5 * PI * t) : 1.0; };
  Eigen::VectorXd mu =
      GaussLobattoParabolic::evolveIBVPGaussLobatto(fe_space, T, M, g);

  // Write the solution to a .vtk file
  const std::string filename = "solution";
  lf::io::VtkWriter vtk_writer(mesh_p, filename + ".vtk");
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  for (int i = 0; i < mu.size(); ++i) {
    nodal_data->operator()(dofh.Entity(i)) = mu(i);
  };
  vtk_writer.WritePointData(filename, *nodal_data);
  std::cout << "Generated " << filename << ".vtk" << std::endl;

  return 0;
}
