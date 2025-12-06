/**
 * @file finitevolumerobin_main.cc
 * @brief NPDE homework FiniteVolumeRobin code
 * @author Philippe Peter
 * @date February 2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "finitevolumerobin.h"

int main() {
  // coefficient functions
  auto g = [](const Eigen::Vector2d& /*x*/) { return 1.0; };
  auto gamma = [](const Eigen::Vector2d& x) { return 1.0 + x(0) * x(0); };

  // The equation is solved on  the four test meshes
  // disk1.msh, disk2.msh, disk3.msh and disk4.msh
  for (int i = 1; i <= 4; ++i) {
    // read mesh
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory),
                              "meshes/disk" + std::to_string(i) + ".msh");
    auto mesh_p = reader.mesh();

    // Construct dofhanlder for linear finite elements on the current mesh.
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

    // Create a dataset of boolean flags indicating edges on the boundary of the
    // mesh
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

    // ASSEMBLE GALERKIN MATRIX
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

    // First the part corresponding to piecewise Lagrangian finite elements
    lf::uscalfe::LinearFELaplaceElementMatrix elmat_provider;
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);

    // Next the part corresponding to the modifications of the Galerkin matrix
    // on the boundary.
    FiniteVolumeRobin::EdgeMatrixProvider edmat_provider(gamma, bd_flags);
    lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edmat_provider, A);

    // RIGHT-HAND SIDE VECTOR
    Eigen::VectorXd phi(dofh.NumDofs());
    phi.setZero();

    // Contributions on the boundary to the rhs vector
    FiniteVolumeRobin::EdgeVectorProvider edvec_provider(g, bd_flags);
    lf::assemble::AssembleVectorLocally(1, dofh, edvec_provider, phi);

    // SOLVE LINEAR SYSTEM
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    Eigen::VectorXd sol_vec = solver.solve(phi);

    // OUTPUT RESULTS TO VTK FILE
    // construct mesh function representing the finite element solution
    lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    // construct vtk writer
    lf::io::VtkWriter vtk_writer(
        mesh_p, "finite_volume_robin_solution_" + std::to_string(i) + ".vtk");
    // output data
    vtk_writer.WritePointData(
        "finite_volume_robin_solution_" + std::to_string(i), mf_sol);
  }
  return 0;
}
