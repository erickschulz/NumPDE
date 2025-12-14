#ifndef STABLE_EVALUATION_AT_A_POINT_H
#define STABLE_EVALUATION_AT_A_POINT_H

/**
 * @file stableevaluationatapoint.h
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher, Erick Schulz & Philippe Peter
 * @date 29.11.2021
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <memory>
#include <utility>

namespace StableEvaluationAtAPoint {

/** @brief Approximates the mesh size for the given mesh.*/
double MeshSize(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p);

/** @brief Returns the outer normal of the unit squre at point x*/
Eigen::Vector2d OuterNormalUnitSquare(Eigen::Vector2d x);

class FundamentalSolution {
 public:
  FundamentalSolution(Eigen::Vector2d x) : x_{x} {}

  // Computes  G_x(y)
  double operator()(Eigen::Vector2d y);
  // Computes the gradient of  G_x(y)
  Eigen::Vector2d grad(Eigen::Vector2d y);

 private:
  Eigen::Vector2d x_;
};

/** @brief Evaluates the Integral P_SL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * @warning The supplied mesh object must hold a triangulation of the **unit
 * square**. This functions only works in this particular setting
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
double PSL(std::shared_ptr<const lf::mesh::Mesh> mesh_p, FUNCTOR&& v,
           const Eigen::Vector2d x) {
  double value = 0.0;
  FundamentalSolution G(x);
#if SOLUTION
  // Flag edges on the boundary
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Loop over boundary edges
  for (const lf::mesh::Entity* e : mesh_p->Entities(1)) {
    if (bd_flags_edge(*e)) {
      const lf::geometry::Geometry* geo_ptr = e->Geometry();
      LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");

      // Fetch coordinates of corner points
      const Eigen::Matrix2d corners = lf::geometry::Corners(*geo_ptr);
      // Determine midpoint of edges
      const Eigen::Vector2d midpoint{0.5 * (corners.col(0) + corners.col(1))};

      // Compute and add the edge contribution
      value += v(midpoint) * G(midpoint) * lf::geometry::Volume(*geo_ptr);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return value;
}
/* SAM_LISTING_END_1 */

/** @brief Evaluates the Integral P_DL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 *
 * @warning The supplied mesh object must hold a triangulation of the **unit
 * square**. This functions only works in this particular setting
 *
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
double PDL(std::shared_ptr<const lf::mesh::Mesh> mesh_p, FUNCTOR&& v,
           const Eigen::Vector2d x) {
  double value = 0.0;
  FundamentalSolution G(x);
#if SOLUTION
  // Flag edges on the boundary
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Loop over boundary edges
  for (const lf::mesh::Entity* e : mesh_p->Entities(1)) {
    if (bd_flags_edge(*e)) {
      const lf::geometry::Geometry* geo_ptr = e->Geometry();
      LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");

      // Fetch coordinates of corner points
      Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
      // Determine midpoints of edges
      const Eigen::Vector2d midpoint{0.5 * (corners.col(0) + corners.col(1))};

      // Determine the normal vector n on the unit square.
      Eigen::Vector2d n = OuterNormalUnitSquare(midpoint);

      // Compute and the elemental contribution
      value += v(midpoint) * (G.grad(midpoint)).dot(n) *
               lf::geometry::Volume(*geo_ptr);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return value;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
/** @brief  This function computes u(x) = P_SL(grad u * n) - P_DL(u).
 * For u(x) = log( (x + (1, 0)^T).norm() ) and x = (0.3, 0.4)^T,
 * it computes the difference between the analytical and numerical
 * evaluation of u.
 *
 * @warning The supplied mesh object must hold a triangulation of the **unit
 * square**. This functions only works in this particular setting.
 */
double PointEval(std::shared_ptr<const lf::mesh::Mesh> mesh_p);
/* SAM_LISTING_END_3 */

class Psi {
 public:
  Psi(Eigen::Vector2d center) : center_(center) {}

  // computes Psi_x(y)
  double operator()(Eigen::Vector2d y);
  // computes grad(Psi_x)(y)
  Eigen::Vector2d grad(Eigen::Vector2d y);
  // computes the laplacian of Psi_x at y
  double lapl(Eigen::Vector2d y);

 private:
  Eigen::Vector2d center_ = Eigen::Vector2d(0.5, 0.5);
};

/** @brief Computes Jstar
 * @param fe_space: finite element space defined on a triangular mesh of the
 * square
 * @param uFE: Coefficient vector of the finite element function wrt the
 * fe_space
 * @param x: Evaluation point
 */
/* SAM_LISTING_BEGIN_4 */
double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
             Eigen::VectorXd uFE, const Eigen::Vector2d x);

/* SAM_LISTING_END_4 */

/** @brief Verifies that the assumptions on Psi_x are satisfied and evaluates
 * Jstar
 * @param fe_space: finite element space defined on a triangular mesh of the
 * square
 * @param uFE: Coefficient vector of the finite element function wrt the
 * fe_space
 * @param x: Evaluation point
 */
double StablePointEvaluation(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    Eigen::VectorXd uFE, const Eigen::Vector2d x);

/** @brief Solves the Laplace equation using Dirichlet conditions g */
template <typename FUNCTOR>
Eigen::VectorXd SolveBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    FUNCTOR&& g) {
  Eigen::VectorXd discrete_solution;

  // Extract mesh and Dofhandler
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  auto N_dofs = dofh.NumDofs();

  // Obtain specification for shape functions on edges
  const auto* rsf_edge_p =
      fe_space_p->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Dirichlet data
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};
  // Right-hand side source function f
  lf::mesh::utils::MeshFunctionConstant mf_f{0.0};

  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // Compute Galerkin Matrix
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Compute right-hand side vector
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Impose essential Boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  auto edges_flag_values_Dirichlet{
      lf::fe::InitEssentialConditionFromFunction(*fe_space_p, bd_flags, mf_g)};
  // Eliminate Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&edges_flag_values_Dirichlet](lf::assemble::glb_idx_t gdof_idx) {
        return edges_flag_values_Dirichlet[gdof_idx];
      },
      A, phi);

  // Assembly completed! Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();

  // II : SOLVING  THE LINEAR SYSTEM
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  discrete_solution = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return discrete_solution;
};
/**
 * @brief Evaluates a finite element function at a point specified by its global
 * coordinates
 */
double EvaluateFEFunction(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& uFE, Eigen::Vector2d global, double tol = 10E-10);

/** @brief Returns the result of evaluating u_h(x) directly or by the stable
 * scheme */
template <typename FUNCTOR>
std::pair<double, double> ComparePointEval(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    FUNCTOR&& g, Eigen::Vector2d x) {
  double direct_eval = 0.0;
  double stable_eval = 0.0;
#if SOLUTION
  // Compute FE solution:
  Eigen::VectorXd uFE = SolveBVP(fe_space, g);

  // use the two evaluation methods:
  direct_eval = EvaluateFEFunction(fe_space, uFE, x);
  stable_eval = StablePointEvaluation(fe_space, uFE, x);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return {direct_eval, stable_eval};
}

} /* namespace StableEvaluationAtAPoint */

#endif  // STABLE_EVALUATION_AT_A_POINT_H