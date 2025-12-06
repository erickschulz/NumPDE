/**
 * @file radauthreetimestepping.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz, edited by Oliver Rietmann
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimestepping.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>

namespace RadauThreeTimestepping {

/**
 * @brief Implementation of the right hand side (time dependent) source vector
 * for the parabolic heat equation
 * @param dofh A reference to the DOFHandler
 * @param time The time at which to evaluate the source vector
 * @returns The source vector at time `time`
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler& dofh,
                                    double time) {
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Right-hand side vector has to be set to zero initially
  Eigen::VectorXd phi(N_dofs);
  // Functor for computing the source function at 2d coordinates
  auto f = [time](Eigen::Vector2d x) -> double {
    const double PI = 3.14159265358979323846;
    Eigen::Vector2d v(std::cos(time * PI), std::sin(time * PI));
    return ((x - 0.5 * v).norm() < 0.5) ? 1.0 : 0.0;
  };
  auto mesh_p = dofh.Mesh();  // pointer to current mesh
  phi.setZero();

  /* Assembling right-hand side source vector */
  // Initialize object taking care of local computations on all cells.
  TrapRuleLinFEElemVecProvider<decltype(f)> elvec_builder(f);
  // Computing right hand side vector
  // Invoke assembly on cells (codim == 0 as first agrument)
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  /* Enforce the zero Dirichlet boundary conditions */
  // Obtain an array of boolean flags for the vertices of the mesh: 'true'
  // indicates that the vertex lies on the boundary. This predicate will
  // guarantee that the computations are carried only on the boundary vertices
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

  // Assigning zero to the boundary values of phi
  for (const lf::mesh::Entity* vertex : mesh_p->Entities(2)) {
    if (bd_flags(*vertex)) {
      auto dof_idx = dofh.GlobalDofIndices(*vertex);
      LF_ASSERT_MSG(
          dofh.NumLocalDofs(*vertex) == 1,
          "Too many global indices were returned for a vertex entity!");
      phi(dof_idx[0]) = 0.0;
    }
  }
  return phi;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Heat evolution solver: the solver obtains the
 * discrete evolution operator from the Radau3MOLTimestepper class and
 * repeatedly iterates its applicaiton starting from the initial condition
 * @param dofh The DOFHandler object
 * @param m is total number of steps until final time final_time (double)
 * @param final_time The duration for which to solve the PDE
 * @returns The solution at the final timestep
 */
/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler& dofh,
                                   unsigned int m, double final_time) {
  Eigen::VectorXd discrete_heat_sol(dofh.NumDofs());
  double tau = final_time / m;                          // step size
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());  // dim. of FE space

  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "\n>>> SolveHeatEvolution: m = " << m << ", N = " << N_dofs
            << std::endl;
  /* Setting up the problem information */
  // Precomputing the required data for the Runge-Kutta method
  // Assemble the Runge-Kutta Radau IIA 2-stages method solver (order 3)
  Radau3MOLTimestepper radau_solver(dofh);
  // Starting with the zero initial condition vector
  Eigen::VectorXd discrete_solution_cur =
      radau_solver.discreteEvolutionOperator(0.0, tau,
                                             Eigen::VectorXd::Zero(N_dofs));

  std::cout << "\n>> Iterating the action of discreteEvolutionOperator"
            << std::endl;
  /* Evolving the parabolic heat system */
  // While less elegant, we use a current and next step solution vector in the
  // iteration to stay away from potential harming aliasing effects of putting
  // an Eigen::Vector on both sides of an assignment statement.
  Eigen::VectorXd discrete_solution_next;
  for (int i = 1; i < m; i++) {
    discrete_solution_next = radau_solver.discreteEvolutionOperator(
        i * tau, tau, discrete_solution_cur);
    discrete_solution_cur = discrete_solution_next;
  }
  discrete_heat_sol = discrete_solution_cur;
  return discrete_heat_sol;
}
/* SAM_LISTING_END_6 */

/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity& tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Compute the area of the triangle cell
  const double area = lf::geometry::Volume(*(tria.Geometry()));
  // Assemble the mass element matrix over the cell
  // clang-format off
  elMat << 2.0, 1.0, 1.0,
           1.0, 2.0, 1.0,
           1.0, 1.0, 2.0;
  // clang-format on
  elMat *= area / 12.0;
  return elMat;  // return the local mass element matrix
}

/* Implementing constructor of class Radau3MOLTimestepper */
/* SAM_LISTING_BEGIN_4 */
Radau3MOLTimestepper::Radau3MOLTimestepper(const lf::assemble::DofHandler& dofh)
    : dofh_(dofh) {
  std::cout << "\n>> Constructing SRadau3MOLTimestepper " << std::endl;
  auto mesh_p = dofh.Mesh();  // pointer to current mesh

  // Instantiating Galerkin matrices to be pre-computed
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Matrices in triplet format holding Galerkin matrices, zero initially.
  lf::assemble::COOMatrix<double> A_COO(N_dofs,
                                        N_dofs);  // element matrix Laplace
  lf::assemble::COOMatrix<double> M_COO(N_dofs,
                                        N_dofs);  // element mass matrix

  std::cout << "> Initializing the Galerking local matrices builders"
            << std::endl;
  // Initialize classes containing the information required for the
  // local computations of the Galerkin matrices. Simple implementations of
  // LinFEMassMatrixProvider and TrapRuleLinFEElemVecProvider adapted to this
  // particular problem was written to spare some of the overhead calculations
  // involved in the use of the more general LehrFEM++ matrices providers.
  lf::uscalfe::LinearFELaplaceElementMatrix elLapMat_builder;
  LinFEMassMatrixProvider elMassMat_builder;

  std::cout << "> Assembling Galerking matrices in COO format" << std::endl;
  // Compute the Galerkin matrices
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrices A and M.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elLapMat_builder, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMassMat_builder, M_COO);

  // Enforcing zero Dirichlet boundary conditions
  // Obtain an array of boolean flags for the vertices of the mesh: 'true'
  // indicates that the vertex lies on the boundary.
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsColumns
  auto bdy_vertices_selector = [&bd_flags, &dofh](unsigned int idx) -> bool {
    return bd_flags(dofh.Entity(idx));
  };
  dropMatrixRowsColumns(bdy_vertices_selector, A_COO);
  dropMatrixRowsColumns(bdy_vertices_selector, M_COO);

  std::cout << "> Converting triplets to sparse matrices" << std::endl;
  // Creating the private Galerkin stiffness and mass matrices
  A_ = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();

  // Runge-Kutta matrices defining the 2-stage Radau timestepping. In the
  // Butcher tableau, this corresponds to c = (1/3 1)^T (top-left column
  // vector), b^T = (3/4 1/4) (bottom-right row vector), U_11 = 5/12, U_12 =
  // -1/12, U_21 = 3/4, U_22 = 1/4 (top-right block); values are fixed in
  // time
  // clang-format off
    U_ << 5.0/12.0, -1.0/12.0,
              0.75,      0.25;
    c_ << 1.0/3.0, 1.0;
    b_ << 0.75, 0.25;
  // clang-format on
  // Precomputing the kronecker products involved in the implicit linear system
  // for the increments of the RADAU-2 method
  M_Kp_ = Eigen::kroneckerProduct(Eigen::Matrix<double, 2, 2>::Identity(), M);
  A_Kp_ = Eigen::kroneckerProduct(U_, A_);
}
/* SAM_LISTING_END_4 */

/* Implementation of Radau3MOLTimestepper member functions */
// The function discreteEvolutionOperator() returns the discretized evolution
// operator as obtained from the Runge-Kutta Radau IIA 2-stages method using the
// Butcher table as stored in the Radau3MOLTimestepper class
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd Radau3MOLTimestepper::discreteEvolutionOperator(
    double time, double tau, const Eigen::VectorXd& mu) const {
  Eigen::VectorXd discrete_evolution_operator(dofh_.NumDofs());
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh_.NumDofs());
  LF_VERIFY_MSG(N_dofs == mu.size(),
                "Dimension mismatch between the number of degrees of freedom "
                "and the dimension of the argument vector.");

  // Building the linear system for the implicitely defined increments
  // Assembling the right hand side using block initialization
  Eigen::VectorXd linSys_rhs(2 * N_dofs);
  Eigen::VectorXd rhs_subtraction_term = A_ * mu;  // precomputation
  linSys_rhs << rhsVectorheatSource(dofh_, time + c_[0] * tau) -
                    rhs_subtraction_term,
      rhsVectorheatSource(dofh_, time + tau) - rhs_subtraction_term;

  // Implicit Runge-Kutta methods lead to systems of equations that must be
  // solved in order to obtained the increments.
  Eigen::SparseMatrix<double> linSys_mat;
  // Assembling the system right hand side matrix using the (unfortunately
  // officially not supported) Eigen Kronecker product
  linSys_mat = M_Kp_ + tau * A_Kp_;
  LF_VERIFY_MSG(linSys_mat.rows() == linSys_mat.cols(),
                "The linSys_mat Eigen matrix is not squared.");
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(linSys_mat);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::VectorXd k_vec = solver.solve(linSys_rhs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  // Compute action of the discrete evolution operator on argument vec
  discrete_evolution_operator = mu + tau * (b_[0] * k_vec.topRows(N_dofs) +
                                            b_[1] * k_vec.bottomRows(N_dofs));
  return discrete_evolution_operator;
}
/* SAM_LISTING_END_5 */

}  // namespace RadauThreeTimestepping
