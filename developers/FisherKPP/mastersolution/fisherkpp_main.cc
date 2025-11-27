/** @file fisherkpp_main.cc
 *  @brief Homework Problem FisherKPP
 *  @author Amélie Justine Loher
 *  @date 12.05.20
 *  @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// In the interest of not changing the problem text over at
// https://gitlab.math.ethz.ch/ralfh/npdeflipped, I will tell clang-tidy to
// ignore the bugprone-suspicious-include warning. (Manuel Saladin, 2024-05-28)
#include "fisherkpp.cc"  // NOLINT(bugprone-suspicious-include)

using namespace FisherKPP;

// Simulation of human migration
void humanmigration();

void humanmigration() {
  // Obtain mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/earth.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  // Finite Element Space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Dofhandler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Initial Population: located in Eritrea
  Eigen::VectorXd u0(N_dofs);
  u0.setZero();
  u0(277) = 80;

  // Taking the topography into account, we have a non-constant diffusioen
  // Coefficient. Mountain chains impede the dispersal of the population.
  Eigen::Vector2d Himalaya(-280, 29);
  Eigen::Vector2d Alps(-350, 46);
  Eigen::Vector2d Karakoram(-285, 36);
  Eigen::Vector2d Hindukush(-288, 36);
  Eigen::Vector2d RockyMountains(-109, 44);
  Eigen::Vector2d Ural(-300, 60);
  Eigen::Vector2d Andes(-66, -21);
  auto c = [&Himalaya, &Alps, &Karakoram, &Hindukush, &RockyMountains, &Ural,
            &Andes](Eigen::Vector2d x) -> double {
    double diffCoeff = 90.0;

    if ((x - Himalaya).norm() <= 3) {
      diffCoeff = 15.0;
      std::cout << "Take Himalaya into account." << std::endl;
    }

    if ((x - Karakoram).norm() <= 3) {
      diffCoeff = 16.0;
      std::cout << "Take Karakoram into account." << std::endl;
    }

    if ((x - Hindukush).norm() <= 2) {
      diffCoeff = 16.0;
      std::cout << "Take Hindukush into account." << std::endl;
    }

    if ((x - Alps).norm() <= 3) {
      diffCoeff = 15.0;
      std::cout << "Take Alps into account." << std::endl;
    }

    if ((x - Ural).norm() <= 2) {
      diffCoeff = 18.0;
      std::cout << "Take Ural into account." << std::endl;
    }

    if ((x - RockyMountains).norm() <= 2) {
      diffCoeff = 20.0;
      std::cout << "Take Himalaya into account." << std::endl;
    }

    if ((x - Andes).norm() <= 3) {
      diffCoeff = 19.0;
      std::cout << "Take Himalaya into account." << std::endl;
    }

    return diffCoeff;
  };

  // Growth Factor
  double lambda = 4.94;

  /* Strang Splitting Method
   * SDIRK-2 evolution of linear parabolic term
   * Exact Evolution for nonlinear reaction term
   */

  // Total number of timesteps
  unsigned int m = 100;
  double T = 1.;  // the timestepsize tau will equal T/m = 0.01

  std::cout << "You are running the simulation on the globe." << std::endl;

  // First we assemble the carrying capacity maps.
  Eigen::MatrixXd car_cap(N_dofs, 16);

  Eigen::VectorXd K(N_dofs);
  K.setZero();
  K = 0.2 * Eigen::VectorXd::Ones(N_dofs);
  auto c_cap = [](Eigen::Vector2d x) -> double { return 80.0; };
  // NOTE: c = 80, lambda = 0, K does not matter.
  StrangSplit DiffusionCapacity(fe_space, T, m, 0.0, c_cap);

  Eigen::VectorXd k0(N_dofs);

  // t = 200 kya - 150 kya
  k0.setZero();
  k0(277) = 80;
  car_cap.col(0) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(1) = DiffusionCapacity.Evolution(K, car_cap.col(0));
  car_cap.col(2) = DiffusionCapacity.Evolution(K, car_cap.col(1));
  car_cap.col(2) *= 10;

  std::cout << "Carrying Capacity 200kya - 150kya!" << std::endl;

  // t = 150 kya - 130 kya
  std::cout << "Carrying Capacity 150kya - 130kya!" << std::endl;

  // t = 130 kya - 100 kya
  k0.setZero();
  k0(334) = 80;
  car_cap.col(3) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(4) = DiffusionCapacity.Evolution(K, car_cap.col(3));
  car_cap.col(5) = DiffusionCapacity.Evolution(K, car_cap.col(4));
  car_cap.col(5) *= 10;

  std::cout << "Carrying Capacity 130kya - 100kya!" << std::endl;

  // t = 100 kya - 70 kya
  std::cout << "Carrying Capacity 100kya - 70kya!" << std::endl;

  // t = 70 kya - 65 kya
  k0.setZero();
  k0(253) = 80;
  car_cap.col(6) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(7) = DiffusionCapacity.Evolution(K, car_cap.col(6));
  car_cap.col(7) *= 15;

  std::cout << "Carrying Capacity 70kya - 65kya!" << std::endl;

  // t = 65 kya - 50 kya
  k0.setZero();
  k0(1775) = 80;
  k0(222) = 80;
  k0(181) = 80;
  car_cap.col(8) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(9) = DiffusionCapacity.Evolution(K, car_cap.col(8));
  car_cap.col(9) *= 8;

  std::cout << "Carrying Capacity 65kya - 50kya!" << std::endl;

  // t = 50 kya - 45 kya
  k0.setZero();
  k0(400) = 200;
  car_cap.col(10) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(10) *= 2;

  std::cout << "Carrying Capacity 50kya - 45kya!" << std::endl;

  // t = 45 kya - 25 kya
  k0.setZero();
  k0(492) = 80;
  k0(114) = 80;
  car_cap.col(11) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(12) = DiffusionCapacity.Evolution(K, car_cap.col(11));
  car_cap.col(12) *= 8;

  std::cout << "Carrying Capacity 45kya - 25kya!" << std::endl;
  // t = 25 kya - 15 kya
  k0.setZero();
  k0(1342) = 80;
  car_cap.col(13) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(14) = DiffusionCapacity.Evolution(K, car_cap.col(13));
  car_cap.col(14) *= 30;

  std::cout << "Carrying Capacity 25kya - 15kya!" << std::endl;

  // t = 15 kya - 1 kya
  k0.setZero();
  k0(1257) = 80;
  k0(1100) = 80;
  car_cap.col(15) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(15) *= 5;

  std::cout << "Carrying Capacity 15kya - 0kya!" << std::endl;
  std::cout << "Capacity maps are assembled!" << std::endl;

  /* Now we may compute the solution */
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c);

  Eigen::MatrixXd sol(N_dofs, 20);
  Eigen::VectorXd cap(N_dofs);
  cap.setZero();
  // t = 200 kya - 150 kya
  cap = car_cap.col(2);
  sol.col(0) = StrangSplitter.Evolution(cap, u0);

  std::cout << "Solution 200kya - 150kya!" << std::endl;

  // t = 150 kya - 130 kya
  sol.col(1) = StrangSplitter.Evolution(cap, sol.col(0));

  std::cout << "Solution 150kya - 130kya!" << std::endl;

  // t = 130 kya - 100 kya
  cap = cap + car_cap.col(5);
  sol.col(2) = StrangSplitter.Evolution(cap, sol.col(1));

  std::cout << "Solution 130kya - 100kya!" << std::endl;

  // t = 100 kya - 70 kya
  sol.col(3) = StrangSplitter.Evolution(cap, sol.col(2));

  std::cout << "Solution 100kya - 70kya!" << std::endl;

  // t = 70 kya - 65 kya
  cap = cap + car_cap.col(7);
  sol.col(4) = StrangSplitter.Evolution(cap, sol.col(3));
  sol.col(5) = StrangSplitter.Evolution(cap, sol.col(4));

  std::cout << "Solution 70kya - 65kya!" << std::endl;

  // t = 65 kya - 50 kya
  cap = cap + car_cap.col(9);
  sol.col(6) = StrangSplitter.Evolution(cap, sol.col(5));
  sol.col(7) = StrangSplitter.Evolution(cap, sol.col(6));

  std::cout << "Solution 65kya - 50kya!" << std::endl;

  // t = 50 kya - 45 kya
  cap = cap + car_cap.col(10);
  sol.col(8) = StrangSplitter.Evolution(cap, sol.col(7));

  std::cout << "Solution 50kya - 45kya!" << std::endl;

  // t = 45 kya - 25 kya
  cap = cap + car_cap.col(12);
  sol.col(9) = StrangSplitter.Evolution(cap, sol.col(8));
  sol.col(10) = StrangSplitter.Evolution(cap, sol.col(9));
  sol.col(11) = StrangSplitter.Evolution(cap, sol.col(10));
  sol.col(12) = StrangSplitter.Evolution(cap, sol.col(11));

  std::cout << "Solution 45kya - 25kya!" << std::endl;

  // t = 25 kya - 15 kya
  cap = cap + car_cap.col(14);
  sol.col(13) = StrangSplitter.Evolution(cap, sol.col(12));
  sol.col(14) = StrangSplitter.Evolution(cap, sol.col(13));

  std::cout << "Solution 25kya - 15kya!" << std::endl;

  // t = 15 kya - 1 kya
  cap = cap + car_cap.col(15);
  sol.col(15) = StrangSplitter.Evolution(cap, sol.col(14));
  sol.col(16) = StrangSplitter.Evolution(cap, sol.col(15));
  sol.col(17) = StrangSplitter.Evolution(cap, sol.col(16));
  sol.col(18) = StrangSplitter.Evolution(cap, sol.col(17));
  sol.col(19) = StrangSplitter.Evolution(cap, sol.col(18));

  std::cout << "Solution 15kya - 1kya!" << std::endl;
  std::cout << "Note that for sake of simplicity, we did not take into account "
               "the dispersal over sea. "
            << "We shall not expect that the islands (especially Australia) "
               "will be inhabited. "
            << "This would have required non local boundary conditions."
            << std::endl;

  for (int k = 1; k < 21; k++) {
    std::stringstream filename;
    filename << "sol" << k << "_human_migration.vtk";

    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) =
          sol.col(k - 1)[global_idx];
    }

    vtk_writer.WritePointData("sol", *nodal_data);
  }
}

void modelproblem();

/* SAM_LISTING_BEGIN_9 */
void modelproblem() {
  std::cout << "You are running the model problem, i.e. you solve the Fisher "
               "equation on the model domain. "
            << std::endl;
  // Obtain mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/island.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  // Finite Element Space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Dofhandler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Initial Population density
  Eigen::VectorXd u0(N_dofs);
  u0.setZero();
  u0(321) = 0.3;
  u0(567) = 0.3;

  // Diffusion Coefficient
  auto c = [](Eigen::Vector2d x) -> double { return 1.2; };
  double lambda = 2.1;  // Growth Factor
  // Carrying capacity
  Eigen::VectorXd K{0.8 * Eigen::VectorXd::Ones(N_dofs)};
  // Time Steps
  unsigned int m = 100;
  double T = 1.;

  // Compute the solution with method of lines and Strang splitting
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c);
  // Five snapshots
  std::vector<Eigen::VectorXd> sol;
  std::cout << "Computing solution after 100 timesteps..." << std::endl;
  sol.push_back(StrangSplitter.Evolution(K, u0));
  // Uncomment the following calls to StrangSplitter.Evolution in order
  // to solve for a longer evolution.
  /*std::cout << "Computing solution after 200 timesteps..." << std::endl;
  sol.push_back(StrangSplitter.Evolution(K, sol[0]));
  std::cout << "Computing solution after 300 timesteps..." << std::endl;
  sol.push_back(StrangSplitter.Evolution(K, sol[1]));
  std::cout << "Computing solution after 400 timesteps..." << std::endl;
  sol.push_back(StrangSplitter.Evolution(K, sol[2]));
  std::cout << "Computing solution after 500 timesteps..." << std::endl;
  sol.push_back(StrangSplitter.Evolution(K, sol[3]));*/

  // Use VTK-Writer for Visualization of solution.
  std::cout << "Writting solution(s) in VTK format." << std::endl;
#if SOLUTION
  std::cout << std::size(sol) << std::endl;
  for (int k = 0; k < std::size(sol); k++) {
    std::stringstream filename;
    filename << "model_problem_sol" << k + 1 << ".vtk";
    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) = sol[k][global_idx];
    }
    vtk_writer.WritePointData("sol", *nodal_data);
  }
  std::cout
      << "Solution after i*100 timesteps written to model_problem_sol'i'.vtk"
      << std::endl;
#else
  //====================
  // Your code goes here: write .vtk files sol[1-6].vtk
  //====================
#endif
}
/* SAM_LISTING_END_9 */

int main(int argc, char** argv) {
  std::cout << "Usage: " << argv[0] << " [h|m]" << std::endl;
  std::cout << "h: human migration simulation" << std::endl;
  std::cout << "m: model problem simulation" << std::endl;
  std::cout << "\nFinite-element simulation of the Fisher/KPP evolution"
            << std::endl;
  std::cout << "===============================================" << std::endl;

  if (argc < 2) {
    std::cout << "No option selected: Will run both..." << std::endl;
    humanmigration();
    modelproblem();
    return 0;
  }

  std::string selection(argv[1]);
  switch (selection[0]) {
    case 'h': {
      humanmigration();
      break;
    }
    case 'm': {
      modelproblem();
      break;
    }
    default: {
      std::cout << "Unrecognized input: terminating .." << std::endl;
      break;
    }
  }
  return 0;
}
