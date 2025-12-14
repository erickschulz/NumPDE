/**
 * @file
 * @brief Solution of source-free heat equation and computation of H1
 *  	  seminorms on different triangular meshes and refinement levels
 * @author Julien Gacon, Amélie Loher
 * @date   March 2019
 */

#include "unstablebvp.h"
// General includes
#include <array>
#include <fstream>
#include <memory>
#include <string>
// Eigen
#include <Eigen/Core>
#include <Eigen/SparseLU>
// Lehrfempp
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

namespace UnstableBVP {

std::shared_ptr<lf::refinement::MeshHierarchy> createMeshHierarchy(
    const int reflevels, const std::string& mesh_type) {
  // Helper object: mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Decide where the triangular domain should be located in x_2 direction
  // by adding an offset to the x_2 coordinate of the nodes
  double offset = 0;
  if (mesh_type == "top") {
    offset = 1.5;
  } else if (mesh_type == "bottom") {
    offset = -1.5;
  } else {
    // already at 0
  }

  // Define the nodes
  std::array<std::array<double, 2>, 3> node_coord{
      std::array<double, 2>({0.5, -0.5 + offset}),
      std::array<double, 2>({0, 0.5 + offset}),
      std::array<double, 2>({1, 0.5 + offset})};

  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(Eigen::Vector2d({node[0], node[1]}));
  }

  // Initialize triangle
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<lf::base::size_type>({0, 1, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));

  // Get a pointer to the mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // (optional) Print information about the mesh
  // std::cout << "   Mesh info\n" << *mesh_p;

  // Ask LehrFEM++ to create a hierarchy of nested meshes
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);

  return multi_mesh_p;
}

double solveTemperatureDistribution(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // **********************************************************************
  // Stage 0: provide all coefficient functions mainly through lambda
  //          functions and derived MeshFunctions
  // **********************************************************************

  // The boundary condition
  auto bc = [](Eigen::Vector2d x) -> double {
    return x[1] <= 0 ? 1 - x[1] : 0;
  };
  // Wrap into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_bc{bc};

  // We use lowest-order (p.w. linear Lagrangian finite elements), for which
  // LehrFEM++ provides a built-in description according to the paradigm of
  // parametric finite elements.
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // **********************************************************************
  // Stage 1: Assemble finite element Galerkin matrix
  // **********************************************************************

  // Dimension of finite element space`
  const lf::base::size_type N_dofs(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  // Element matrix builder for the negative Laplacian
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};

  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // **********************************************************************
  // Stage 2: Right-hand side vector
  // **********************************************************************

  // Define RHS vector
  // No source, hence it is simply zero
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // **********************************************************************
  // Stage 3: Fixing solution components according to essential (Dirichlet)
  //          boundary conditions
  // **********************************************************************

  // Obtain specification for shape functions on edges
  // NOLINTBEGIN(clang-analyzer-deadcode.DeadStores)
  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
  // NOLINTEND(clang-analyzer-deadcode.DeadStores)
  LF_ASSERT_MSG(rsf_edge_p != nullptr, "FE specification for edges missing");

  // Obtain an array of boolean flags for the edges (codim 1) of the mesh,
  // `true` indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // edges.
  auto ess_bdc_flags_values{lf::fe::InitEssentialConditionFromFunction(
      *fe_space,
      [&bd_flags](const lf::mesh::Entity& edge) -> bool {
        return (bd_flags(edge));
      },
      mf_bc)};

  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::assemble::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  // **********************************************************************
  // Stage 4: Solve LSE
  // **********************************************************************

  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  // **********************************************************************
  // Stage 5: Compute H1 seminorm
  // **********************************************************************

  // construct a mesh function that represent |grad(x)|^2 and integrate it over
  // the mesh to get the H^1 semi-norm:
  auto mf_norm = squaredNorm(lf::fe::MeshFunctionGradFE(fe_space, sol_vec));
  const double norm = std::sqrt(IntegrateMeshFunction(*mesh_p, mf_norm, 2));

  return norm;
}

}  // namespace UnstableBVP
