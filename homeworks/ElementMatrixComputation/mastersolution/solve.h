/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 06.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef SOLVE_H_
#define SOLVE_H_

#include <lf/assemble/assemble.h>
#include <lf/assemble/assemble_concepts.h>
#include <lf/fe/fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "../meshes/mesh.h"
#include "mylinearfeelementmatrix.h"
#include "mylinearloadvector.h"

namespace ElementMatrixComputation {

/**
 * @brief      Given element providers, it creates FE space, assembles Galerkin
 * matrix and vector, and solves the LSE.
 *
 * @param[in]  mesh_p         A mesh pointer
 * @param[in]  elmat_provider  An element matrix provider object
 * @param[in]  elvec_provider  An element vector provider object
 *
 * @return     The solution vector
 */
/* SAM_LISTING_BEGIN_1 */
template <lf::assemble::EntityMatrixProvider ELMAT_PROVIDER,
          lf::assemble::EntityMatrixProvider ELVEC_PROVIDER>
Eigen::VectorXd solve(ELMAT_PROVIDER& elmat_provider,
                      ELVEC_PROVIDER& elvec_provider) {
  // Use one of LehrFEM++'s default meshes. Try different meshes by changing the
  // function index parameter. See the documentation of that function for
  // details ablut the available meshes
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(8, 1.0 / 3.0);
  // We use a linear Lagrangian FE space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Reference to current mesh, obtained from the FE space
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space`
  const lf::base::size_type N_dofs(dofh.NumDofs());

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Invoke assembly on cells (co-dimension = 0). The element matrix provider is
  // passed as an argument
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Invoke assembly on cells (codim == 0). The element vector provider is
  // passed as an argument
  AssembleVectorLocally(0, dofh, elvec_provider, phi);

  // Define solution vector
  Eigen::VectorXd sol_vec = Eigen::VectorXd::Zero(N_dofs);

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  // Make sure LU-decomposition could be computed
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  // Backward substitution
  sol_vec = solver.solve(phi);
  /* SAM_LISTING_END_1 */

  double solver_error = (A_crs * sol_vec - phi).norm();
  double solver_relative_error = solver_error / phi.norm();
  std::cout << "Error of solver: Absolute error = " << solver_error
            << ", relative error = " << solver_relative_error << std::endl;

  // Postprocessing: Compute and output norms of the finite element solution
  double L2norm = std::sqrt(IntegrateMeshFunction(
      *mesh_p, squaredNorm(lf::fe::MeshFunctionFE(fe_space, sol_vec)), 2));
  double H1snorm = std::sqrt(IntegrateMeshFunction(
      *mesh_p, squaredNorm(lf::fe::MeshFunctionGradFE(fe_space, sol_vec)), 2));

  std::cout << "Norms of FE solution: L2-norm = " << L2norm
            << ", H1-seminorm = " << H1snorm << std::endl;
  return sol_vec;
}

/**
 * @brief      Right hand side functional `f(x)`
 * @param[in]  x     A 2D vector that provides the evaluation of the source
 * @return     The value of `f` evaluated at point `x`
 */
inline double f(Eigen::Vector2d x) { return 1 + x(0) * x(0) + x(1) * x(1); };

/**
 * @brief Solve Poisson's equation -△u = f using LehrFEM++'s built in element
 * matrix and vector providers
 * @return     The solution vector
 */
Eigen::VectorXd solvePoissonBVP();

/**
 * @brief Solve Neumann equation -△u + u = f where ɑ is a constant diffusion
 * coefficient using a custom implementation of element matrix and element
 * vector providers
 * @return     The solution vector
 */
Eigen::VectorXd solveNeumannEq();

}  // namespace ElementMatrixComputation

#endif  // SOLVE_H_
