/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "sdirkmethodoflines.h"

namespace SDIRKMethodOfLines {

/** @brief This class implements a Lehrfem++ matrix provider defining a
LinFEMassMatrixProvider::Eval function returning the local MASS matrix for
linear first-order lagrange FE bases over triangular mesh (only!). Integration
over the triangular cells is performed using the trapezoidal rule.*/
class LinFEMassMatrixProvider {
 public:
  /** @brief default constructor */
  explicit LinFEMassMatrixProvider() = default;
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell refers to current cell for which the element vector is desired
   * The implementation uses an analytic formula defined over triangular cells*/
  Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity& tria);
};  // class LinFEMassMatrixProvider
/** Implementing member function Eval of class LinFEMassMatrixProvider*/
/* SAM_LISTING_BEGIN_8 */
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity& tria) {
  Eigen::Matrix<double, 3, 3> elMat;
#if SOLUTION
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Compute the area of the triangle cell
  const double area = lf::geometry::Volume(*(tria.Geometry()));
  // Compute the mass element matrix over the cell
  // clang-format off
  elMat << 2.0, 1.0, 1.0,
           1.0, 2.0, 1.0,
           1.0, 1.0, 2.0;
  // clang-format on
  elMat *= area / 12.0;
#else
  //====================
  // Your code goes here
  //====================
#endif

  return elMat;  // return the local mass element matrix
}  // LinFEMassMatrixProvider::Eval
/* SAM_LISTING_END_8 */

/** @brief This class implements a Lehrfem++ matrix provider defining a
LinearMassEdgeMatrixProvider<FUNCTOR>::Eval function returning the local EDGE
MASS matrix for linear first-order lagrange FE over triangular mesh (only!).
Integration over the triangular cells is performed using the trapezoidal rule.*/
template <typename FUNCTOR>
class LinearMassEdgeMatrixProvider {
 public:
  /** @brief Constructor storing the right hand side function
   *  @param predicate is lambda function returning booleans on edge entities
   *         cool_coeff is the convective cooling coefficient */
  explicit LinearMassEdgeMatrixProvider(FUNCTOR predicate, double cool_coeff)
      : predicate_(predicate), cool_coeff_(cool_coeff) {}
  /** @brief Default implement: all edges are active */
  virtual bool isActive(const lf::mesh::Entity& /*edge*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param edge is current entity for which the element vector is desired
   * The implementation uses simple vertex based quadrature */
  Eigen::Matrix<double, 2, 2> Eval(const lf::mesh::Entity& edge);

 private:
  /** predicate_ provides booleans for boundary edges */
  FUNCTOR predicate_;
  double cool_coeff_;
};  // class LinearMassEdgeMatrixProvider
/* Implementing member function Eval of class LinearMassEdgeMatrixProvider */
/* SAM_LISTING_BEGIN_9 */
template <typename FUNCTOR>
Eigen::Matrix<double, 2, 2> LinearMassEdgeMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& edge) {
  Eigen::Matrix<double, 2, 2> elBdyEdgeMat;
#if SOLUTION
  // Throw error in case not edge entity
  LF_VERIFY_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                "Unsupported entity type " << edge.RefEl());
  if (predicate_(edge)) {
    // Obtain endpoints of the edge
    auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
    LF_ASSERT_MSG(endpoints.cols() == 2, "Wrong no endpoints in " << edge);
    // Compute length of the edge
    double edge_length = (endpoints.col(1) - endpoints.col(0)).norm();
    // Compute the local edge mass element matrix
    // clang-format off
    elBdyEdgeMat << 1.0 / 3.0, 1.0 / 6.0, 
                    1.0 / 6.0, 1.0 / 3.0;
    // clang-format on
    elBdyEdgeMat *= edge_length * cool_coeff_;

    return elBdyEdgeMat;  // return the local edge mass element matrix
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return (Eigen::Matrix<double, 2, 2>::Zero());
}  // LinearMassEdgeMatrixProvider<FUNCTOR>::Eval
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler& dofh,
                         double cool_coeff) {
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      sparse_pair;
#if SOLUTION
  std::cout << ">> Constructing SDIRK2Timestepper " << std::endl;
  auto mesh_p = dofh.Mesh();  // pointer to current mesh
  // Instantiating Galerkin matrices to be pre-computed
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Obtain an array of boolean flags for the edges of the mesh: 'true'
  // indicates that the edge lies on the boundary. This predicate will
  // guarantee that the computations are carried only on the boundary edges
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Creating predicate that will guarantee that the computations are carried
  // only on the edges of the mesh using the boundary flags
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge);
  };
  // Matrices in triplet format holding Galerkin matrices, zero initially.
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);

  std::cout << "> Initializing the Galerking local matrices builders"
            << std::endl;
  // Initialize classes containing the information required for the local
  // computations of the Galerkin matrices. Simple implementations of
  // LinFEMassMatrixProvider adapted to this particular problem was written to
  // spare some of the overhead calculations involved in the use of the more
  // general LehrFEM++ matrices providers.
  lf::uscalfe::LinearFELaplaceElementMatrix elLapMat_builder;
  LinearMassEdgeMatrixProvider<decltype(edges_predicate)> el_MassEdge_builder(
      edges_predicate, cool_coeff);
  LinFEMassMatrixProvider elMassMat_builder;

  std::cout << "> Assembling Galerking matrices in COO format" << std::endl;
  // Compute the Galerkin matrices
  // Invoke assembly on cells and edges (co-dimension = 0 and 1 respectively as
  // first argument). Information about the mesh and the local-to-global map is
  // passed through a Dofhandler object, argument 'dofh'. This function call
  // adds triplets to the internal COO-format representation of sparse matrices
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elLapMat_builder, A_COO);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, el_MassEdge_builder,
                                      A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMassMat_builder, M_COO);

  std::cout << "> Converting to triplets to sparse matrices" << std::endl;
  // Assembling the private Galerkin matrices
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();
  sparse_pair = std::make_pair(A, M);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return sparse_pair;
}
/* SAM_LISTING_END_1 */

/* Implementation of class SDIRK2Timestepper */
// Implementation of SDIRK2Timestepper constructor
#if SOLUTION
SDIRK2Timestepper::SDIRK2Timestepper(const lf::assemble::DofHandler& dofh,
                                     double tau /*nb. steps*/,
                                     double cool_coeff /*cooling coeff*/)
    : tau_(tau), lambda_(1.0 - 0.5 * sqrt(2.0)) {
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      sparse_pair = assembleGalerkinMatrices(dofh, cool_coeff);
  A_ = sparse_pair.first;
  Eigen::SparseMatrix<double> M = sparse_pair.second;

  std::cout << "> Precomputing LU decompositions for Eigen solvers"
            << std::endl;
  // Implicit Runge-Kutta methods lead systems of equations that must be solved
  // in order to obtained the increments. For time evolution problem with fixed
  // step size, the matrix associated to these linear systems is independent of
  // time. We precompute the sparse LU solvers for efficiency.

  solver.compute(M + tau * lambda_ * A_);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
#else
SDIRK2Timestepper::SDIRK2Timestepper(const lf::assemble::DofHandler& dofh,
                                     double tau /*nb. steps*/,
                                     double cool_coeff /*cooling coeff*/)
    : tau_(tau) {
//====================
// Your code goes here
//====================
#endif
}  // SDIRK2Timestepper constructor

/* Implementation of SDIRK2Timestepper member function */
/* SAM_LISTING_BEGIN_9 */
Eigen::VectorXd SDIRK2Timestepper::discreteEvolutionOperator(
    const Eigen::VectorXd& mu) const {
  Eigen::VectorXd discrete_evolution_operator;
#if SOLUTION
  Eigen::VectorXd rhs_vec = -A_ * mu;  // precomputation
  // Stage 1 of SDIRK-2
  Eigen::VectorXd k1_vec = solver.solve(rhs_vec);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // Stage 2 of SDIRK-2
  Eigen::VectorXd k2_vec =
      solver.solve(rhs_vec - tau_ * (1 - lambda_) * A_ * k1_vec);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");

  discrete_evolution_operator =
      mu + tau_ * (1 - lambda_) * k1_vec + tau_ * lambda_ * k2_vec;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return discrete_evolution_operator;
}  // SDIRK2Timestepper::discreteEvolutionOperator
/* SAM_LISTING_END_9 */

/** @Brief Implementing the temperature evolution solver. The solver obtains the
discrete evolution operator from the SDIRK2Timestepper class and repeatedly
iterates its applicaiton starting from the initial condition argument
* @param cool_coeff is the convective cooling coefficient */
/* SAM_LISTING_BEGIN_6 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> solveTemperatureEvolution(
    const lf::assemble::DofHandler& dofh, unsigned int m, double cool_coeff,
    Eigen::VectorXd initial_temperature_vec) {
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solution_pair;
#if SOLUTION
  double tau = 1.0 / m;                                 // step size
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());  // dim. of FE space
  std::cout << "\n>>> SolveTemperatureEvolution: m = " << m
            << ", N = " << N_dofs << ", cool_coeff= " << cool_coeff << "\n"
            << std::endl;

  // Initialization for the energies of the solution at times t = 0, tau, 2*tau,
  // etc... This information is not required to solve the parabolic evolution
  Eigen::VectorXd energies(m + 1);

  // Precomputing the required data for SDIRK-2
  SDIRK2Timestepper SDIRK2_stepper(dofh, tau, cool_coeff);

  std::cout << "\n>> Iterating the action of discreteEvolutionOperator"
            << std::endl;
  // Starting the evolution using the initial conditions
  Eigen::VectorXd discrete_solution_cur =
      SDIRK2_stepper.discreteEvolutionOperator(initial_temperature_vec);
  energies[0] = thermalEnergy(dofh, initial_temperature_vec);
  energies[1] = thermalEnergy(dofh, discrete_solution_cur);
  // Evolving the parabolic temperature system
  // While less elegant, we use a current and next step solution vector in
  // the iteration to stay away from potential harming aliasing effects of
  // putting an Eigen::Vector on both sides of an assignment statement.
  Eigen::VectorXd discrete_solution_next;
  for (int i = 1; i < m; i++) {
    discrete_solution_next =
        SDIRK2_stepper.discreteEvolutionOperator(discrete_solution_cur);
    discrete_solution_cur = discrete_solution_next;
    energies[i + 1] = thermalEnergy(dofh, discrete_solution_cur);
  }
  solution_pair = std::make_pair(discrete_solution_cur, energies);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return solution_pair;
}  // solveTemperatureEvolution

/* SAM_LISTING_END_6 */

/** @brief thermalEnergy computes the integral of the discrete argument function
passed as a vector of coefficients whose index is the global index of the
degree of freedom (linear lagrange basis function) that they multiply */
/* SAM_LISTING_BEGIN_7 */
double thermalEnergy(const lf::assemble::DofHandler& dofh,
                     const Eigen::VectorXd& temperature_vec) {
  double thermal_energy = 0.0;
#if SOLUTION
  auto mesh_p = dofh.Mesh();  // pointer to mesh

  // The thermal energy of the system is defined as the integral of the
  // temperature solution over the domain. It is computed for linear lagrangian
  // finite elements using the trapezoidal rule by summing up the contribution
  // of that quadrature rule over each triangle
  double thermal_energy_loc;
  for (const lf::mesh::Entity* tria : mesh_p->Entities(0)) {
    thermal_energy_loc = 0.0;
    // Compute the area of the triangle
    const double area = lf::geometry::Volume(*(tria->Geometry()));
    // Obtain the global indices of the nodal degrees of freedom
    for (const lf::assemble::gdof_idx_t& g_idx : dofh.GlobalDofIndices(*tria)) {
      thermal_energy_loc += temperature_vec[g_idx];
    }
    // Sum local contribution of quadrature rule
    thermal_energy += thermal_energy_loc * area / 3.0;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return thermal_energy;
}  // thermalEnergy

/* SAM_LISTING_END_7 */

}  // namespace SDIRKMethodOfLines
