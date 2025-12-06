/**
 * @file mixedfemwave.h
 * @brief NPDE homework MixedFEMWave
 * @author Erick Schulz
 * @date 14.06.2020
 * @copyright Developed at ETH Zurich
 */

#ifndef MIXEDFEMWAVE_H_
#define MIXEDFEMWAVE_H_

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>

namespace MixedFEMWave {

// An auxiliary function providing a QuadRule object for the local trapezoidal
// rule on a triangle
lf::quad::QuadRule make_TriaQR_TrapezoidalRule();

/* SAM_LISTING_BEGIN_1 */
template <typename FFUNCTION>
Eigen::VectorXd computeRHS(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V,
    FFUNCTION&& f, double t) {
  // TOOLS AND DATA
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh_V{fe_space_V->LocGlobMap()};
  // Pointer to mesh
  auto mesh_p = dofh_V.Mesh();
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());
  // Right-hand side source function f
  auto mf_f = lf::mesh::utils::MeshFunctionGlobal(
      [f, t](Eigen::Vector2d x) -> double { return f(x, t); });

  // ASSEMBLY
  // Right hand side vector, initialized with zeros
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhs_vec(N_dofs_V);
  rhs_vec.setZero();
  // Computing right-hand side vector using cell-oriented assembly based on the
  // local trapezoidal rule. The information about the quadrature rule to be
  // used can be passed to the constructor through an associative array, which
  // is built first
  // Lehrfempp requires the specification of a quadrature rule for
  // quadrilateral. A generic rule for QUAD that is NOT used in this problem is
  // passed. Alternatively : see computeRHS_alt

  // The following initialization will be enough for the next release of
  // LehrFEM++ lf::uscalfe::quad_rule_collection_t qr_collection{
  //    {lf::base::RefEl::kTria(), make_TriaQR_TrapezoidalRule()}};
  // So far we need a dummy quadrature rule for quadrilaterals
  lf::uscalfe::quad_rule_collection_t qr_collection{
      {lf::base::RefEl::kTria(), make_TriaQR_TrapezoidalRule()},
      {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_P4O4()}};
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_V, mf_f, qr_collection);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh_V, elvec_builder, rhs_vec);

  // Node-index array of flags marking boundary nodes
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Set entries of the RHS vector to zero
  for (int i = 0; i < N_dofs_V; i++) {
    if (bd_flags(dofh_V.Entity(i))) {
      rhs_vec(i) = 0;
    }
  }

  return rhs_vec;
}
/* SAM_LISTING_END_1 */

// Auxiliary function: Determine combined areas of cells adjacent to the nodes
// of a mesh
lf::mesh::utils::CodimMeshDataSet<double> areasOfAdjacentCells(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

// Another way to implement the function: used for debugging
/* SAM_LISTING_BEGIN_4 */
template <typename FFUNCTION>
Eigen::VectorXd computeRHS_alt(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V,
    FFUNCTION&& f, double t) {
  // The current mesh object
  const lf::mesh::Mesh& mesh{*(fe_space_V->Mesh())};
  // Flag nodes on the boundary
  const auto bd_node_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(fe_space_V->Mesh(), 2)};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh_V{fe_space_V->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());
  // Right hand side vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhs_vec(N_dofs_V);
  // Determine area of mesh neighborhood of every node
  const auto areas{areasOfAdjacentCells(fe_space_V->Mesh())};
  // Loop over all nodes
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    // Find global index of dof associated with the node
    const std::span<const lf::assemble::gdof_idx_t> node_dofs(
        dofh_V.GlobalDofIndices(*node));
    LF_ASSERT_MSG(node_dofs.size() == 1, "A node has to carry exatly one dof!");
    // Global index of tent function sitting on the node
    const lf::assemble::gdof_idx_t dof_no = node_dofs[0];
    // If the node is located on the boundary, set entry of right-hand-side
    // vector to zero
    if (bd_node_flags(*node)) {
      rhs_vec[dof_no] = 0.0;
    } else {  // Interior node
      // Obtain position vector of node
      const Eigen::Vector2d node_pos{
          lf::geometry::Corners(*(node->Geometry())).col(0)};
      // Evaluate source function at node
      const double f_val = f(node_pos, t);
      // Directly set entry of right-hand-side vector
      rhs_vec[dof_no] = areas(*node) * f_val / 3.0;
    }
  }
  return rhs_vec;
}
/* SAM_LISTING_END_4 */

Eigen::SparseMatrix<double> computeMQ(const lf::assemble::DofHandler& dofh_Q);

/* SAM_LISTING_BEGIN_6 */
template <typename RHOFUNCTION>
Eigen::SparseMatrix<double> computeMV(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V,
    RHOFUNCTION&& rho) {
  // TOOLS AND DATA
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_V->Mesh();
  const lf::assemble::DofHandler& dofh_V{fe_space_V->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());
  // For returning the matrix
  Eigen::SparseMatrix<double> M_V(N_dofs_V, N_dofs_V);
  auto zero_mf = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  auto rho_mf = lf::mesh::utils::MeshFunctionGlobal(
      [rho](Eigen::Vector2d x) -> double { return rho(x); });

  // ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> M_V_COO(N_dofs_V, N_dofs_V);
  // ReactionDiffusionElementMatrixProvider<SCALAR,DIFF_COEFF,REACTION_COEFF,
  // quadrule>
  // Tell this ENTITY_MATRIX_PROVIDER object to use the local trapezoidal rule
  // Lehrfempp requires the specification of a quadrature rule for
  // quadrilateral. A generic rule for QUAD that is NOT used in this problem is
  // passed. Alternatively : see computeMV_alt
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(zero_mf),
                                                      decltype(rho_mf)>
      elmat_builder(fe_space_V, zero_mf, rho_mf,
                    {{lf::base::RefEl::kTria(), make_TriaQR_TrapezoidalRule()},
                     {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_P4O4()}});
  lf::assemble::AssembleMatrixLocally(0, dofh_V, dofh_V, elmat_builder,
                                      M_V_COO);
  M_V = M_V_COO.makeSparse();
  std::cout << "Assembly: M_V finished" << std::endl;

  return M_V;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
template <typename RHOFUNCTION>
Eigen::SparseMatrix<double> computeMV_alt(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V,
    RHOFUNCTION&& rho) {
  const lf::mesh::Mesh& mesh{*(fe_space_V->Mesh())};
  const lf::assemble::DofHandler& dofh_V{fe_space_V->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());

  // Reserve space for a sparse diagonal matrix
  Eigen::SparseMatrix<double> M_V(N_dofs_V, N_dofs_V);
  M_V.reserve(Eigen::VectorXi::Constant(N_dofs_V, 1));

  // Determine area of mesh neighborhood of every node
  const auto areas{areasOfAdjacentCells(fe_space_V->Mesh())};
  // Loop over all nodes
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    // Obtain position vector of node
    const Eigen::Vector2d node_pos{
        lf::geometry::Corners(*(node->Geometry())).col(0)};
    // Evaluate coefficient function at node
    const double rho_val = rho(node_pos);
    // Find global index of dof associated with the node
    const std::span<const lf::assemble::gdof_idx_t> node_dofs(
        dofh_V.GlobalDofIndices(*node));
    LF_ASSERT_MSG(node_dofs.size() == 1, "A node has to carry exatly one dof!");
    // Global index of tent function sitting on the node
    const lf::assemble::gdof_idx_t dof_no = node_dofs[0];
    M_V.insert(dof_no, dof_no) = areas(*node) * rho_val / 3.0;
  }
  return M_V;
}
/* SAM_LISTING_END_7 */

class BElemMatProvider {
 public:
  explicit BElemMatProvider() = default;
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  Eigen::Matrix<double, 2, 3> Eval(const lf::mesh::Entity& tria);
};

Eigen::SparseMatrix<double> computeB(const lf::assemble::DofHandler& dofh_V,
                                     const lf::assemble::DofHandler& dofh_Q);

/* SAM_LISTING_BEGIN_L */
template <typename RHOFUNCTION, typename FFUNCTION,
          typename RECORDER = std::function<void(double, double)>>
std::pair<Eigen::VectorXd, Eigen::VectorXd> leapfrogMixedWave(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_V,
    const lf::assemble::UniformFEDofHandler& dofh_Q, RHOFUNCTION&& rho,
    FFUNCTION&& f, double T, unsigned int nb_timesteps,
    RECORDER&& rec = [](double, double) {}) {
  // Size of timestep
  double stepsize = T / nb_timesteps;
  // Index mapper for linear Lagrangian FE space
  const lf::assemble::DofHandler& dofh_V{fe_space_V->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs_Q(dofh_Q.NumDofs());
  const lf::uscalfe::size_type N_dofs_V(dofh_V.NumDofs());

  // Zero INITIAL CONDITIONS
  Eigen::VectorXd mu_init = Eigen::VectorXd::Zero(N_dofs_V);
  Eigen::VectorXd kappa_init = Eigen::VectorXd::Zero(N_dofs_Q);

// PRECOMPUTATIONS: essential for efficiency
  std::cout << "Precomputing Galerkin Matrices :" << std::endl;
  // Galerkin matrix $\wt{\VM}_V$
  Eigen::SparseMatrix<double> M_V = computeMV(fe_space_V, rho);
  std::cout << "Creating SparseLU solver for M_V: ";
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_MV;
  solver_MV.compute(M_V);
  LF_VERIFY_MSG(solver_MV.info() == Eigen::Success, "LU decomposition failed");
  std::cout << " Done" << std::endl;
  // Galerkin matrix $\wt{\VM}_Q$
  Eigen::SparseMatrix<double> M_Q = computeMQ(dofh_Q);
  std::cout << "Creating SparseLU solver for M_Q: ";
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_MQ;
  solver_MQ.compute(M_Q);
  LF_VERIFY_MSG(solver_MQ.info() == Eigen::Success, "LU decomposition failed");
  std::cout << "Done" << std::endl;
  // Galerkin matrix $\wt{\VB}$
  Eigen::SparseMatrix<double> B = computeB(dofh_V, dofh_Q);

  // EVOLUTION
  // Vectors for returning result
  Eigen::VectorXd mu;
  Eigen::VectorXd kappa;
  std::cout << "Evolution using the leapfrog timestepping scheme." << std::endl;
  Eigen::VectorXd mu_cur = mu_init;        // $\vec{\mubf}^{(j)}$
  Eigen::VectorXd kappa_cur = kappa_init;  // $\vec{\kappabf}^{(j+\frac12)}$
  Eigen::VectorXd mu_next;                 // $\vec{\mubf}^{(j+1)}$
  Eigen::VectorXd kappa_next;              // $\vec{\kappabf}^{(j+\frac32)}$
  Eigen::VectorXd kappa_avg;
  for (int j = 0; j < nb_timesteps; j++) {
    // Right hand side vector $\vec{\varphibf}(\tau(j+\frac12))$
    const Eigen::VectorXd rhs{computeRHS(fe_space_V, f, stepsize * (j + 0.5))};
    // Update of $\vec{\mubf}$
    mu_next = solver_MV.solve(stepsize * (rhs + B.transpose() * kappa_cur) +
                              M_V * mu_cur);
    // Update of $\vec{\kappabf}$
    kappa_next = solver_MQ.solve(M_Q * kappa_cur - stepsize * B * mu_next);
    // Compute approximate total energy according to \prbeqref{eq:en}
    // "Kinetic energy" at time $\tau(j+1)$
    const double en1 = 0.5 * mu_next.dot(M_V * mu_next);
    // "Potential energy" at time $\tau(j+1)$
    kappa_avg = 0.5 * (kappa_cur + kappa_next);
    const double en2 = 0.5 * kappa_avg.dot(M_Q * kappa_avg);
    // Record contributions to energy, see \prbautoref{sp:a}
    rec(en1, en2);
    // Prepare next timestep
    mu_cur = mu_next;
    kappa_cur = kappa_next;
  }  // end timestepping loop
  mu = mu_cur;
  kappa = kappa_avg;
  return {mu, kappa};
}
/* SAM_LISTING_END_L */

}  // namespace MixedFEMWave

#endif
