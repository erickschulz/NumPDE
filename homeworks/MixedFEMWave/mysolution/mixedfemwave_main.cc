/**
 * @file mixedfemwave_main.cc
 * @brief NPDE homework MixedFEMWave
 * @author Erick Schulz
 * @date 14.06.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <math.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include "mixedfemwave.h"
#include "systemcall.h"

using namespace MixedFEMWave;

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "\n ****** Problem MixedFEMWave ******" << std::endl;
  // PROBLEM DATA
  double T = 1.5;
  unsigned int nb_timesteps = T * 500;
  auto f = [](Eigen::Vector2d x, double t) -> double {
    if (std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 2.5, 2) < 0.25) {
      return 15.0 * ((t < 0.5) ? std::sin(2.0 * lf::base::kPi * t) : 0.0);
    }
    return 0.0;
  };
  auto rho = [](Eigen::Vector2d x) -> double { return 1.0; };

  // LOADING MESH
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 "meshes/bassin3.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  // FINITE ELEMENT SPACES AND DOFHS
  // Scalar finite element space for lowest-order Lagrangian finite elements
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Scalar dofhandler as built along with the finite-element space for V
  const lf::assemble::DofHandler& dofh_V = fe_space_V->LocGlobMap();
  // Dimension of unconstrained finite-element space V
  lf::base::size_type N_dofs_V = dofh_V.NumDofs();
  // Vector dofhandler for the finite element space Q
  lf::assemble::UniformFEDofHandler dofh_Q(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 0},
                                            {lf::base::RefEl::kSegment(), 0},
                                            {lf::base::RefEl::kTria(), 2},
                                            {lf::base::RefEl::kQuad(), 2}});
  // Array for recording energies
  Eigen::VectorXd energy(nb_timesteps);
  // Recorder lambda function: risky, because we rely on the fact that it is not
  // called more than nb_timesteps times
  int cnt = 0;
  auto recorder = [&energy, &cnt](double en1, double en2) {
    energy[cnt++] = en1 + en2;
  };

  // Evolving the wave equation
  auto [sol_u, sol_j] =
      leapfrogMixedWave(fe_space_V, dofh_Q, rho, f, T, nb_timesteps, recorder);

  // Output energy results to csv file
  const Eigen::VectorXd t{Eigen::VectorXd::LinSpaced(nb_timesteps + 1, 0.0, T)};
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream energy_csv;
  energy_csv.open("energy.csv");
  energy_csv << t.tail(nb_timesteps).transpose().format(CSVFormat) << std::endl;
  energy_csv << energy.transpose().format(CSVFormat) << std::endl;
  energy_csv.close();
  std::cout << "Generated energy.csv" << std::endl;
  systemcall::execute("python3 scripts/plot_energy.py energy.csv energy.eps");

  // Output wave solution results to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "wave_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs_V; global_idx++) {
    nodal_data->operator()(dofh_V.Entity(global_idx)) = sol_u[global_idx];
  };
  vtk_writer.WritePointData("wave_solution", *nodal_data);
}
