/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "outputimpedancebvp.h"

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

namespace OutputImpedanceBVP {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveImpedanceBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::Vector2d g) {
  // Related implementations:
  // Homework problem ErrorEstimatesForTraces:
  // https://gitlab.math.ethz.ch/ralfh/npdecodes/tree/master/homeworks/ErrorEstimatesForTraces

  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  Eigen::VectorXd discrete_solution(N_dofs);

  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // I.i : Computing volume matrix for negative Laplace operator
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_2 */
  // I.ii : Computing mass edge matrix resulting from Robin B.C.
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Creating a predicate that will guarantee that the computations are carried
  // only on the interior boundary edges of the mesh using the boundary flags
  auto edges_predicate_RobinBC =
      [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    if (bd_flags(edge)) {
      auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
      return endpoints(0, 0) > 0.05 && 0.95 > endpoints(0, 0) &&
             endpoints(1, 0) > 0.05 && 0.95 > endpoints(1, 0);
    }
    return false;
  };
  // Coefficients used in the class template
  // MassEdgeMatrixProvider< SCALAR, COEFF, EDGESELECTOR >
  auto eta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(eta),
                                      decltype(edges_predicate_RobinBC)>
      edgemat_builder(fe_space_p, eta, edges_predicate_RobinBC);
  // Invoke assembly on edges by specifying co-dimension = 1
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);
  /* SAM_LISTING_END_2 */

  /* SAM_LISTING_BEGIN_9 */
  // I.iii : Computing right-hand side vector
  // Right-hand side source function f
  auto mf_f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // I.iv : Imposing essential boundary conditions
  // Dirichlet data
  auto mf_g = lf::mesh::utils::MeshFunctionGlobal(
      [&g](Eigen::Vector2d x) -> double { return g.dot(x); });
  // Inspired by the example in the documentation of
  // InitEssentialConditionFromFunction()
  // https://craffael.github.io/lehrfempp/namespacelf_1_1uscalfe.html#a5afbd94919f0382cf3fb200c452797ac
  // Creating a predicate that will guarantee that the computations are carried
  // only on the exterior boundary edges of the mesh using the boundary flags
  auto edges_predicate_Dirichlet =
      [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    if (bd_flags(edge)) {
      auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
      if (endpoints(0, 0) <= 0.05 || 0.95 <= endpoints(0, 0) ||
          endpoints(1, 0) <= 0.05 || 0.95 <= endpoints(1, 0)) {
        return true;
      }
    }
    return false;
  };
  // Determine the fixed dofs on the boundary and their values
  // Alternative: See lecturedemoDirichlet() in
  // https://github.com/craffael/lehrfempp/blob/master/examples/lecturedemos/lecturedemoassemble.cc
  auto edges_flag_values_Dirichlet{lf::fe::InitEssentialConditionFromFunction(
      *fe_space_p, edges_predicate_Dirichlet, mf_g)};
  // Eliminate Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&edges_flag_values_Dirichlet](lf::assemble::glb_idx_t gdof_idx) {
        return edges_flag_values_Dirichlet[gdof_idx];
      },
      A, phi);

  // Assembly completed! Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();

// II : SOLVING  THE LINEAR SYSTEM
  // II.i : Setting up Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // II.ii : Solving
  discrete_solution = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

// do nothing
  return discrete_solution;
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_3 */
double computeBoundaryOutputFunctional(
    const Eigen::VectorXd eta,
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::Vector2d d) {
  double func_val = 0.0;
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};

  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Creating a predicate that will guarantee that the computations are carried
  // only on the interior boundary edges of the mesh using the boundary flags
  auto edges_predicate_RobinBC =
      [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    if (bd_flags(edge)) {
      auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
      return endpoints(0, 0) > 0.05 && 0.95 > endpoints(0, 0) &&
             endpoints(1, 0) > 0.05 && 0.95 > endpoints(1, 0);
    }
    return false;
  };

  // Computing value of the functional
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    if (edges_predicate_RobinBC(*edge)) {
      // Find the endpoints global indices
      auto dof_idx = dofh.GlobalDofIndices(*edge);
      assert(dofh.NumLocalDofs(*edge) == 2);
      // Value of linear function in he endpoints
      const double nu0 = eta(dof_idx[0]);
      const double nu1 = eta(dof_idx[1]);
      // Find coordinates of endpoints
      auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
      // Length of edge
      const double elen = (endpoints.col(1) - endpoints.col(0)).norm();
      LF_ASSERT_MSG(endpoints.cols() == 2, "Wrong no endpoints in " << *edge);
      func_val += elen / 6.0 *
                  d.dot(2 * endpoints.col(0) * nu0 + endpoints.col(0) * nu1 +
                        endpoints.col(1) * nu0 + 2 * endpoints.col(1) * nu1);
    }
  }
  return func_val;
};
/* SAM_LISTING_END_3 */

}  // namespace OutputImpedanceBVP
