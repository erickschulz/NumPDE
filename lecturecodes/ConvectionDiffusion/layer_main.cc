/**
 * @file layer_main.cc
 * @brief Solves CD BVP for exact solution with an internal layer
 * @author Philippe Peter
 * @date July 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <memory>
#include <string>

#include "cd_tools.h"
#include "standard_fem.h"
#include "supg.h"
#include "systemcall.h"
#include "upwind.h"

int main() {
  // parameter functions:
  // boundary conditions
  const auto g = [](const Eigen::Vector2d& x) {
    return x(0) > x(1) ? 1.0 : 0.0;
  };
  // velocity field
  const auto v = [](const Eigen::Vector2d& x) {
    return Eigen::Vector2d(1.0, 1.0);
  };
  // diffusion coefficient
  const auto eps = [](const Eigen::Vector2d& x) { return 10E-10; };
  // source function
  const auto f = [](const Eigen::Vector2d& x) { return 0.0; };

  // Read Mesh from file
  std::string mesh_file = "meshes/mesh_square.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();

  // Construct dofhanlder for linear finite element space on the mesh.
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solutions using Standard FE, Upwind, SUPG method
  Eigen::VectorXd sol_standard =
      ConvectionDiffusion::SolveCDBVPStandardFem(fe_space, eps, v, f, g);
  lf::fe::MeshFunctionFE sol_standard_mf(fe_space, sol_standard);

  Eigen::VectorXd sol_stable =
      ConvectionDiffusion::SolveCDBVPUpwind(fe_space, eps, v, f, g);
  lf::fe::MeshFunctionFE sol_upwind_mf(fe_space, sol_stable);

  Eigen::VectorXd sol_supg =
      ConvectionDiffusion::SolveCDBVPSupg(fe_space, eps, v, f, g);
  lf::fe::MeshFunctionFE sol_supg_mf(fe_space, sol_supg);

  // Output solution along the curve gamma
  auto gamma = [](double t) { return Eigen::Vector2d(t, 1 - t); };
  ConvectionDiffusion::SampleMeshFunction("results_standard_FEM.txt", mesh_p,
                                          gamma, sol_standard_mf, 300);
  ConvectionDiffusion::SampleMeshFunction("results_upwind.txt", mesh_p, gamma,
                                          sol_upwind_mf, 300);
  ConvectionDiffusion::SampleMeshFunction("results_supg.txt", mesh_p, gamma,
                                          sol_supg_mf, 300);

  // Plot
  systemcall::execute("python3 scripts/plot_layer.py .");
  return 0;
}