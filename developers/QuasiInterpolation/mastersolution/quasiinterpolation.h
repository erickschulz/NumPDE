/**
 * @file quasiinterpolation.h
 * @brief NPDE exam TEMPLATE HEADER FILE
 * @author Oliver Rietmann
 * @date 15.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef QUASIINTERPOLATION_H_
#define QUASIINTERPOLATION_H_

#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace QuasiInterpolation {

/**
 * @brief Computes the length of the longest edge
 *
 * @param edges range of edges
 * @return length of a longest edge in the range and zero if the range is empty
 */
double maxLength(const std::span<const lf::mesh::Entity* const>& edges);

/**
 * @brief Produces a mapping of vertices $p$ to the corresponding $(K_p,j)$,
 * where $K_p$ is a largest triangle containing $p$ and $j\in\{0,1,2\}$ is the
 * local index of $p$ in $K_p$.
 *
 * @param mesh_p shared pointer to a triangular mesh
 * @return data set storing a mapping as described above
 */
lf::mesh::utils::CodimMeshDataSet<
    std::pair<const lf::mesh::Entity*, unsigned int>>
findKp(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/**
 * @brief Compute the quasi-interpolation operator $Q_h(v)$.
 *
 * @param fe_space FEM space of tent functions to project onto
 * @param v_mf mesh function representing $v(x)$, such that MESHFUNCTION
 * satisfies the LehrFEM concept MeshFunction of value type R=\texttt{double}
 * @return coefficient vector of size equal to the number of DOFs in fe_space,
 * representing the coefficients of the projected function $Q_h(v)$
 */
/* SAM_LISTING_BEGIN_1 */
template <typename MESHFUNCTION>
Eigen::VectorXd quasiInterpolate(
    const lf::uscalfe::FeSpaceLagrangeO1<double>& fe_space,
    MESHFUNCTION&& v_mf) {
  // Get quadrature points and weights on the reference triangle for a
  // quadrature rule that integrates quadratic polynomials exactly
  lf::quad::QuadRule quadrule =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2);
  const int P = quadrule.NumPoints();
  Eigen::MatrixXd zeta_hat = quadrule.Points();
  Eigen::VectorXd w_hat = quadrule.Weights();
#if SOLUTION
  // Precomputation: Cache values of $(\hat{\psi}^1,\hat{\psi}^2,\hat{\psi}^3)$
  // at reference quadrature nodes. Use the formulas \prbeqref{eq:psigen}.
  Eigen::Matrix<double, 3, 2> A;
  A << -24.0, -24.0, 24.0, 0.0, 0.0, 24.0;
  Eigen::Vector3d b{18.0, -6.0, -6.0};
  auto psi_hat = [A, b](Eigen::Vector2d x) -> Eigen::Vector3d {
    return A * x + b;
  };
  Eigen::MatrixXd psi_hat_values(3, P);
  for (int l = 0; l < P; ++l) {
    psi_hat_values.col(l) = psi_hat(zeta_hat.col(l));
  }
#else
  // ========================================
  // First part of your solution here
  // ========================================
#endif
  // Retrieve the map $p \mapsto (K_p, j)$
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space.Mesh();
  lf::mesh::utils::CodimMeshDataSet<
      std::pair<const lf::mesh::Entity*, unsigned int>>
      Kp_mesh_data_set = findKp(mesh_p);

  const lf::assemble::DofHandler& dofh{fe_space.LocGlobMap()};
  const int N = dofh.NumDofs();
  // Vector for returning the basis expansion coefficients of the
  // quasi-interpolant
  Eigen::VectorXd coefficients(N);
#if SOLUTION
  // Compute projection onto each basis function
  for (lf::assemble::gdof_idx_t i = 0; i < N; ++i) {
    // Get $p$, $K_p$ and the local index $j$ of $p$ in $K_p$
    const lf::mesh::Entity& point = dofh.Entity(i);
    std::pair<const lf::mesh::Entity*, unsigned int> Kp_localindex =
        Kp_mesh_data_set(point);
    const lf::mesh::Entity* Kp = Kp_localindex.first;
    unsigned int j = Kp_localindex.second;
    // Evaluate v on quadrature points in $K_p$
    const std::vector<double> v_zeta = v_mf(*Kp, zeta_hat);
    // Compute integral on current triangle Kp
    double sum = 0.0;
    for (int l = 0; l < P; ++l) {
      sum += w_hat(l) * v_zeta[l] * psi_hat_values(j, l);
    }
    coefficients(i) = sum;
  }
#else
  // ========================================
  // Second part of your solution here
  // ========================================
#endif
  return coefficients;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute the interpolation error for convergence plots given a sequence
 * of $L$ refinements.
 *
 * @param l2_error vector of size $L$ to hold the errors w.r.t. $\lVert\ \cdot\
 * \rVert_{L^2(\Omega)}$.
 * @param h1_error vector of size $L$ to hold the errors w.r.t. $\lVert\ \cdot\
 * \rVert_{H^1(\Omega)}$.
 * @param mesh_hierarchy_p sequence of $L$ refinements of a triangular mesh
 * @param u_mf mesh function representing the function $u(x)$ to project,
 * such that MESHFUNCTION_U satisfies the LehrFEM concept MESHFUNCTION of value
 * type R=\texttt{double}
 * @param grad_u_mf mesh function representing the gradient $\nabla u(x)$ of
 * $u(x)$, such that MESHFUNCTION_GRAD_U satisfies the LehrFEM concept
 * MESHFUNCTION of value type R=\texttt{Eigen::Vector2d}
 */
template <typename MESHFUNCTION_U, typename MESHFUNCTION_GRAD_U>
std::pair<Eigen::VectorXd, Eigen::VectorXd> interpolationError(
    std::shared_ptr<lf::refinement::MeshHierarchy> mesh_hierarchy_p,
    MESHFUNCTION_U&& mf_u, MESHFUNCTION_GRAD_U&& mf_grad_u) {
  // Make sure that the coefficient types are compatible
  // static_assert(lf::mesh::utils::MeshFunction<MESHFUNCTION_U>);
  // static_assert(lf::mesh::utils::MeshFunction<MESHFUNCTION_GRAD_U>);
  // Allocate sequences for storing error norms
  int L = mesh_hierarchy_p->NumLevels();
  Eigen::VectorXd l2_error(L);
  Eigen::VectorXd h1_error(L);
  // Main loop over all levels
  for (int k = 0; k < L; ++k) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = mesh_hierarchy_p->getMesh(k);
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute basis expansion coefficient vector of
    Eigen::VectorXd coefficients = quasiInterpolate(*fe_space, mf_u);
    // Compute error norms as in ellbvp_linfe_demo.cc
    // Mesh function representing the solution
    const lf::fe::MeshFunctionFE mf_intp(fe_space, coefficients);
    // Mesh function encoding the gradient of the finite element function
    const lf::fe::MeshFunctionGradFE mf_grad_intp(fe_space, coefficients);
    // Compute the L2 norm of the interpolation error
    l2_error(k) = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_intp - mf_u), 4));
    double h1semi_error = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_intp - mf_grad_u), 4));

    auto square = [](double x) { return x * x; };
    h1_error(k) = std::sqrt(square(h1semi_error) + square(l2_error(k)));
  }
  return {l2_error, h1_error};
}

}  // namespace QuasiInterpolation

#endif  // #ifndef QUASIINTERPOLATION_H_
