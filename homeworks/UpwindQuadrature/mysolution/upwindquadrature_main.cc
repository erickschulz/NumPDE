/**
 * @file upwindquadrature_main.cc
 * @brief NPDE homework template main
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>

#include "../../../lecturecodes/ConvectionDiffusion/convection_emp.h"
#include "upwindquadrature.h"

int main() {
  // PARAMETERS
  // mesh specification (number of cells in both sides of the tensor-product
  // triangular mesh)
  int M = 49;

  // coefficient functions:
  // Dirichlet functor
  const auto g = [](const Eigen::Vector2d& x) {
    return x(1) == 0 ? 0.5 - std::abs(x(0) - 0.5) : 0.0;
  };
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};

  // velocity field
  const auto v = [](const Eigen::Vector2d& x) {
    return (Eigen::Vector2d() << -x(1), x(0)).finished();
  };

  // diffusion coefficient
  const double eps = 1e-4;
  lf::mesh::utils::MeshFunctionConstant mf_eps{eps};

  // MESH CONSTRUCTION
  // construct a triangular tensor product mesh on the unit square
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(M)
      .setNumYCells(M);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

  // DOF HANDLER & FINITE ELEMENT SPACE
  // Construct dofhanlder for linear finite elements on the mesh.
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // PREPARING DATA TO IMPOSE DIRICHLET CONDITIONS
  // Create a dataset of boolean flags indicating edges on the boundary of the
  // mesh
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // boundary.
  auto ess_bdc_flags_values{
      lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g)};

  //============================================================================
  // SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (STANDARD: UNSTABLE)
  //============================================================================
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

  // ASSEMBLE GALERKIN MATRIX
  // First the part corresponding to the laplacian
  lf::uscalfe::ReactionDiffusionElementMatrixProvider laplacian_provider(
      fe_space, mf_eps, lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider, A);

  // Next part corresponding to the convection term:
  ConvectionDiffusion::ConvectionElementMatrixProvider convection_provider(v);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider, A);

  // RIGHT-HAND SIDE VECTOR
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();

  // IMPOSE DIRICHLET CONDITIONS:
  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::uscalfe::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  // SOLVE LINEAR SYSTEM
  Eigen::SparseMatrix A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  // OUTPUT RESULTS TO VTK FILe
  // construct mesh function representing finite element solution
  lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
  // construct vtk writer
  lf::io::VtkWriter vtk_writer(mesh_p,
                               "upwind_quadrature_solution_unstable.vtk");
  // output data
  vtk_writer.WritePointData("upwind_quadrature_solution_unstable", mf_sol);

//============================================================================
// SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (UPWIND: STABLE)
//============================================================================
/* SAM_LISTING_BEGIN_7 */
  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_7 */
  return 0;
}
