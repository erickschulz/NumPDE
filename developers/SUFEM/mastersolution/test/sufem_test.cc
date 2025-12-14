/**
 * @file sufem_test.cc
 * @brief NPDE homework SUFEM code
 * @author R. Hiptmair
 * @date July 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../sufem.h"

#include <gtest/gtest.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <cstddef>

namespace SUFEM::test {

TEST(sufem, matsum) {
  std::cout << "Test: zero row sums of advection Galerkin matrix" << std::endl;

  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs = dofh.NumDofs();
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  SUFEM::AdvectionElementMatrixProvider elmat_provider(
      fe_space, lf::mesh::utils::MeshFunctionGlobal(
                    [](Eigen::Vector2d x) -> Eigen::Vector2d {
                      return Eigen::Vector2d(-x[1], x[0]);
                    }));
  // Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  // Invoke cell-oriented assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A_COO);
  //  Assembly completed: Convert COO matrix A into CRS format
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  // Check row sum
  Eigen::MatrixXd A_dense = A;
  std::cout << "Galerkin matrix = " << A_dense << std::endl;
  EXPECT_NEAR(A_dense.rowwise().sum().norm(), 0.0, 1.0E-10)
      << "Row sum not zero";
}

TEST(sufem, eval) {
  std::cout << "Test: values of advection Galerkin matrix" << std::endl;

  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs = dofh.NumDofs();
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  SUFEM::AdvectionElementMatrixProvider elmat_provider(
      fe_space, lf::mesh::utils::MeshFunctionGlobal(
                    [](Eigen::Vector2d x) -> Eigen::Vector2d {
                      return Eigen::Vector2d(-x[1], x[0]);
                    }));
  // Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);

  // Invoke cell-oriented assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A_COO);
  //  Assembly completed: Convert COO matrix A into CRS format
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  // Check row sum
  Eigen::MatrixXd A_dense = A;

  // We check some random values, note that hardcoding will be checked for and
  // is not allowed
  EXPECT_NEAR(A_dense(0, 1), -0.0162037, 1e-4);
  EXPECT_NEAR(A_dense(10, 12), 0., 1e-4);
  EXPECT_NEAR(A_dense(8, 3), -0.037037, 1e-4);
  EXPECT_NEAR(A_dense(12, 11), 0.116898, 1e-4);
}

TEST(sufem, MeshFunctionDiffTensor) {
  // Velocity field: rigid body rotation
  auto velo = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  // We need the velocity in MF format
  lf::mesh::utils::MeshFunctionGlobal mf_velo(velo);

  MeshFunctionDiffTensor mf_diff(mf_velo);
  std::cout << "Test: values of MeshFunctionDiffTensor" << std::endl;

  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // Initialize correct solutions: note that hardcoding is not allowed and will
  // be checked for!
  std::vector<int> indices = {5, 3, 8};
  std::vector<Eigen::Matrix2d> sols(3);
  sols.at(0) << 0.022222, -0.022222, -0.022222, 0.022222;
  sols.at(1) << 0., 0., 0., 0.0790569;
  sols.at(2) << 0.0707107, -0.0942809, -0.0942809, 0.125708;

  // Compare the output to the correct solution
  for (int i = 0; i < indices.size(); ++i) {
    auto e = mesh_p->Entities(0)[indices.at(i)];
    EXPECT_LE((mf_diff(*e, Eigen::Vector2d(0., 0.)).at(0) - sols.at(i)).norm(),
              1e-6);
  }
}

TEST(sufem, sumat) {
  std::cout << "Test: zero row sums of SU Galerkin matrix" << std::endl;

  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // Set up global FE space; lowest order Lagrangian finite elements
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto A_COO = SUFEM::buildSUGalerkinMatrix(
      fe_space, lf::mesh::utils::MeshFunctionGlobal(
                    [](Eigen::Vector2d x) -> Eigen::Vector2d {
                      return Eigen::Vector2d(-x[1], x[0]);
                    }));
  // Convert COO matrix A into CRS format
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  // Check row sum
  Eigen::MatrixXd A_dense = A;
  std::cout << "Galerkin matrix = " << A_dense << std::endl;
  EXPECT_NEAR(A_dense.rowwise().sum().norm(), 0.0, 1.0E-10)
      << "Row sum not zero";
}

TEST(sufem, inflow) {
  std::cout << "Test, if right nodes on inflow boundary are flagged"
            << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto inflow_flags = SUFEM::flagNodesOnInflowBoundary(
      mesh_p, lf::mesh::utils::MeshFunctionGlobal(
                  [](Eigen::Vector2d x) -> Eigen::Vector2d {
                    return Eigen::Vector2d(-x[1], x[0]);
                  }));
  // Run through all nodes of the mesh, print their location and the flag
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    // Fetch location of the node
    Eigen::Vector2d pos = lf::geometry::Corners(*(node->Geometry())).col(0);
    std::cout << "Node @ [" << pos.transpose()
              << "]: " << (inflow_flags(*node) ? "INFLOW" : "neutral")
              << std::endl;
  }
}

TEST(sufem, convergencerate) {
  unsigned int reflevels = 6;
  const char* filename = "rotation";
  std::cout << "Convergence test for SU FEM on unit square, pure advection"
            << std::endl;
  // Velocity field: rigid body rotation
  auto velo = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  // Exact solution
  auto u = [](Eigen::Vector2d x) -> double {
    const double xn = x.norm();
    return ((xn <= 1.0) ? (1.0 - xn) * xn : 0.0);
  };
  // Initialize mesh function
  lf::mesh::utils::MeshFunctionGlobal mf_velo(velo);
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);

  // Triangular mesh of the unit square
  auto cmesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Obtain a pointer to a hierarchy of nested meshes
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(cmesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes used in test\n";
  multi_mesh.PrintInfo(std::cout);
  // Number of levels
  lf::base::size_type L = multi_mesh.NumLevels();
  // LEVEL LOOP: Do computations on all levels
  // Vector for storing error norms
  std::vector<std::pair<lf::base::size_type, double>> errs{};
  for (lf::base::size_type level = 0; level < L; ++level) {
    std::shared_ptr<lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs = dofh.NumDofs();
    // Solve Dirichlet boundary value problem
    Eigen::VectorXd uh_coeffs =
        solveAdvectionDirichlet(fe_space, mf_velo, mf_u);
    // Mesh function for FE solution
    const lf::fe::MeshFunctionFE mf_uh(fe_space, uh_coeffs);
    // Compute error
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_uh - mf_u), 2));
    errs.emplace_back(N_dofs, L2err);
    // Output result to file
    if (filename) {
      lf::io::VtkWriter vtk_writer(
          mesh_p, std::string(filename) + std::to_string(level) + ".vtk");
      auto nodal_data =
          lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
      for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
        nodal_data->operator()(dofh.Entity(global_idx)) = uh_coeffs[global_idx];
      };
      vtk_writer.WritePointData("solution_advection_problem", *nodal_data);
    }
  }
  std::cout << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 error" << std::endl;
  for (int i = 1; i < errs.size(); ++i) {
    auto [N, l2err] = errs[i];
    auto [Nprev, l2errprev] = errs[i - 1];
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::endl;

    // We expect a convergence rate of 1
    EXPECT_LE(std::log(l2err / l2errprev) / log(2.), -1.);
  }
}

}  // namespace SUFEM::test
