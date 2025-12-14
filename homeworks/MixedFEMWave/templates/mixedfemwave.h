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
  // =============================================
  // Your code here
  // =============================================

  return M_V;
}
/* SAM_LISTING_END_6 */


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
// ========================================
// Your code here
// ========================================

  // EVOLUTION
  // Vectors for returning result
  Eigen::VectorXd mu;
  Eigen::VectorXd kappa;
  // ========================================
  // Your code here
  // ========================================
  return {mu, kappa};
}
/* SAM_LISTING_END_L */

}  // namespace MixedFEMWave

#endif
