/**
 * @file upwindfinitevolume_main.cc
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/gmsh_reader.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "upwindfinitevolume.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  // Read in mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), "meshes/mesh.msh");

  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 6)};
  int num_meshes = mesh_seq_p->NumLevels();

  double eps = 1e-6;

  auto v = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(1, 0.5 - x[1]);
  };

  auto u = [](Eigen::Vector2d x) -> double {
    return std::sin(4 * M_PI * x[1]) * std::exp(-x[0]);
  };

  auto f = [&u, eps](Eigen::Vector2d x) -> double {
    return eps * (1.0 - 16.0 * M_PI * M_PI) * u(x) - 2 * u(x) +
           (0.5 - x[1]) * 4 * M_PI * std::cos(4.0 * M_PI * x[1]) *
               std::exp(-x[0]);
  };

  // Convergence study
  std::vector<double> vec_l2error;
  std::vector<double> vec_Ndofs;
  for (int level = 0; level < num_meshes; ++level) {
    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);

    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(
        cur_mesh, {{lf::base::RefEl::kPoint(), 1},
                   {lf::base::RefEl::kSegment(), 0},
                   {lf::base::RefEl::kTria(), 0},
                   {lf::base::RefEl::kQuad(), 0}});
    int N_dofs = cur_dofh.NumDofs();

    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    UpwindFiniteVolume::ElementMatrixProvider my_mat_provider(v, eps);
    lf::assemble::AssembleMatrixLocally(0, cur_dofh, cur_dofh, my_mat_provider,
                                        A);

    Eigen::VectorXd phi(N_dofs);
    phi.setZero();
    UpwindFiniteVolume::ElementVectorProvider my_vec_provider(f);
    lf::assemble::AssembleVectorLocally(0, cur_dofh, my_vec_provider, phi);

    // Dirichlet boundary condition
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(cur_mesh, 2)};
    auto my_selector = [&cur_dofh, &u, &bd_flags](unsigned int dof_idx) {
      if (bd_flags(cur_dofh.Entity(dof_idx))) {
        const lf::mesh::Entity& entity{cur_dofh.Entity(dof_idx)};
        const lf::geometry::Geometry* geo_p = entity.Geometry();
        const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);
        Eigen::Vector2d corner = corners.col(0);
        return std::make_pair(true, u(corner));
      }
      return std::make_pair(false, 42.0);
    };
    lf::assemble::FixFlaggedSolutionComponents<double>(my_selector, A, phi);

    const Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // Numerical Solution
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    if (solver.info() != Eigen::Success) {
      printf("Decomposition Failed.\n");
    }
    Eigen::VectorXd sol_vec = solver.solve(phi);

    // Exact Solution
    Eigen::VectorXd exact_sol(N_dofs);
    for (int i = 0; i < N_dofs; ++i) {
      const lf::mesh::Entity& entity{cur_dofh.Entity(i)};
      const lf::geometry::Geometry* geo_p = entity.Geometry();
      const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);
      exact_sol[i] = u(corners.col(0));
    }

    // L2 Error
    double l2_error = (sol_vec - exact_sol).norm() / N_dofs;
    vec_l2error.push_back(l2_error);
    vec_Ndofs.push_back(N_dofs);
  }

  for (int i = 0; i < vec_l2error.size(); ++i) {
    if (i > 0) {
      double conv_rate =
          (std::log(vec_l2error.at(i - 1)) - std::log(vec_l2error.at(i))) /
          (std::log(vec_Ndofs.at(i)) - std::log(vec_Ndofs.at(i - 1)));
      std::cout << "N_dofs: " << vec_Ndofs.at(i)
                << "; L2 Error: " << vec_l2error.at(i)
                << "; Conv. rate: " << conv_rate << std::endl;
    } else {
      std::cout << "N_dofs: " << vec_Ndofs.at(i)
                << "; L2 Error: " << vec_l2error.at(i) << std::endl;
    }
  }

  /* SAM_LISTING_END_1 */
  return 0;
}
