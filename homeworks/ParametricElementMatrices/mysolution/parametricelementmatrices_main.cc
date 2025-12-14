/** @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans, Erick Schulz (refactoring)
 * @date 13/03/2019, 19/11/2019 (refactoring)
 * @copyright Developed at ETH Zurich */

// Lehrfempp includes
#include <lf/fe/fe.h>
#include <lf/io/io.h>
// ParametricElementMatrices internal includes
#include "anisotropicdiffusionelementmatrixprovider.h"
#include "fesourceelemvecprovider.h"
#include "impedanceboundaryedgematrixprovider.h"

int main() {
  /* PRODUCE MESH AND FE SPACE */
  // Create a Lehrfem++ square tensor product mesh using internal routines
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(50)
      .setNumYCells(50);
  auto mesh_p = builder.Build();

  // Finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::base::size_type N_dofs(dofh.NumDofs());

  /* GENERATE BVP DATA */
  // Interpolate variable coefficient function w(x) = sin(|x|)
  auto w_func = [](Eigen::Vector2d x) -> double { return std::sin(x.norm()); };
  lf::mesh::utils::MeshFunctionGlobal mf_w_func{w_func};
  auto w = lf::fe::NodalProjection<double>(*fe_space, mf_w_func);
  // Create direction vector entering the anisotropic diffusion tensor
  auto d = [](Eigen::Vector2d x) -> Eigen::Vector2d { return x; };

  /* PRODUCE LINEAR SYSTEM OF EQUATIONS */
  // Initialize Galerkin matrix in triplets format
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Initialize load vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Create the classes containing the information required for the local
  // computations of the Galerkin matrices and the load vector.
  auto elmat_builder =
      ParametricElementMatrices::AnisotropicDiffusionElementMatrixProvider(d);
  auto edgemat_builder =
      ParametricElementMatrices::ImpedanceBoundaryEdgeMatrixProvider(fe_space,
                                                                     w);
  auto elvec_builder =
      ParametricElementMatrices::FESourceElemVecProvider(fe_space, w);
  // Compute the Galerkin matrices
  // Invoke assembly on cells (co-dimension = 0 as first argument) and edges
  // (co-dimension = 1 as first argument). Information about the mesh and the
  // local-to-global map is passed through the Dofhandler object, argument
  // 'dofh'. This function call adds triplets to the internal COO-format
  // representation of the sparse matrixA.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);
  // Compute the load vector
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  /* SOLVE LINEAR SYSTEM OF EQUATIONS */
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  /* SAVE RESULTS */
  // Output results to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p,
                               "parametric_element_matrices_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) = sol_vec[global_idx];
  };
  vtk_writer.WritePointData("parametric_element_matrices_solution",
                            *nodal_data);
  /* SAM_LISTING_END_1 */
  std::cout << "\n The approximate solution was written to:" << std::endl;
  std::cout << ">> parametric_element_matrices_solution.vtk\n" << std::endl;

  return 0;
}
