/**
 * @file robotnavigation.cc
 * @brief NPDE homework RobotNavigation code
 * @author Erick Schulz
 * @date 28.03.2026
 * @copyright Developed at ETH Zurich
 */

#include "robotnavigation.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/entity.h>

namespace RobotNavigation {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solvePoissonBVP(const lf::io::GmshReader& reader) {
#if SOLUTION
  auto mesh_p = reader.mesh();

  // Set up linear Lagrange FE space
  using FESpace = lf::uscalfe::FeSpaceLagrangeO1<double>;
  auto fe_space = std::make_shared<FESpace>(mesh_p);
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  const auto N_dofs = dofh.NumDofs();

  // Assemble Galerkin stiffness matrix (Laplacian)
  constexpr int cell_codim = 0;
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix elmat;
  lf::assemble::AssembleMatrixLocally(cell_codim, dofh, dofh, elmat, A_COO);

  // Assemble RHS load vector with f = 1
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
  auto f = lf::mesh::utils::MeshFunctionConstant(1.0);
  lf::uscalfe::ScalarLoadElementVectorProvider elvec(fe_space, f);
  lf::assemble::AssembleVectorLocally(cell_codim, dofh, elvec, phi);

  constexpr int edge_codim = 1;
  const auto door_id = reader.PhysicalEntityName2Nr("doors", edge_codim);
  auto is_door_edge = [&reader, id = door_id](const lf::mesh::Entity& edge) {
    return reader.IsPhysicalEntity(edge, id);
  };

  // Flag all nodes that lie on a door edge (Dirichlet boundary)
  constexpr int node_codim = 2;
  using NodeFlags = lf::mesh::utils::CodimMeshDataSet<bool>;
  NodeFlags is_door_node{mesh_p, node_codim, false};
  for (auto* edge : mesh_p->Entities(edge_codim)) {
    if (is_door_edge(*edge)) {
      constexpr int endpoint_subcodim = 1;
      const auto& endpoints = edge->SubEntities(endpoint_subcodim);
      for (const auto* node : endpoints) {
        is_door_node(*node) = true;
      }
    }
  }

  // Enforce zero Dirichlet boundary conditions on door nodes
  auto dirichlet_bc =
      [&](lf::assemble::gdof_idx_t dof_idx) -> std::pair<bool, double> {
    const lf::mesh::Entity& node = dofh.Entity(dof_idx);
    return {is_door_node(node), 0.0};
  };
  lf::assemble::FixFlaggedSolutionComponents<double>(dirichlet_bc, A_COO, phi);

  // Solve the linear system
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  return solver.solve(phi);
#else
  return Eigen::VectorXd::Zero(1);
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const lf::mesh::Entity& cell) {
#if SOLUTION
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  auto vertices = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
#else
  return Eigen::Matrix<double, 2, 3>::Zero();
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd GradientProjectionMassMatrixProvider::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
#if SOLUTION
  // 6x6 mass matrix for V_h on a triangle
  Eigen::MatrixXd elMat = Eigen::MatrixXd::Zero(6, 6);
  const double area = lf::geometry::Volume(*(cell.Geometry()));
  // Local basis {(b1, 0), (0, b1), ...,(b3, 0), (0, b3)}
  // clang-format off
  elMat << 2.0, 0.0, 1.0, 0.0, 1.0, 0.0,
           0.0, 2.0, 0.0, 1.0, 0.0, 1.0,
           1.0, 0.0, 2.0, 0.0, 1.0, 0.0,
           0.0, 1.0, 0.0, 2.0, 0.0, 1.0,
           1.0, 0.0, 1.0, 0.0, 2.0, 0.0,
           0.0, 1.0, 0.0, 1.0, 0.0, 2.0;
  // clang-format on
  elMat *= area / 12.0;
  return elMat;
#else
  return Eigen::MatrixXd::Zero(6, 6);
#endif
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd GradientProjectionRhsVectorProvider::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
#if SOLUTION
  const lf::assemble::DofHandler& dofh = fe_space_p_->LocGlobMap();
  auto gdof_idx = dofh.GlobalDofIndices(cell);

  // Gradient of barycentric coordinates
  Eigen::Matrix<double, 2, 3> grad_bary = gradbarycoordinates(cell);

  // Piecewise-constant gradient of u_h on this cell
  Eigen::Vector2d gradu_h = Eigen::Vector2d::Zero();
  for (int i = 0; i < 3; i++) {
    gradu_h += grad_bary.col(i) * u_h_(gdof_idx[i]);
  }

  // Element vector using trapezoidal rule
  const double area = lf::geometry::Volume(*(cell.Geometry()));
  Eigen::VectorXd elVec(6);
  // clang-format off
  elVec << gradu_h(0),
           gradu_h(1),
           gradu_h(0),
           gradu_h(1),
           gradu_h(0),
           gradu_h(1);
  // clang-format on
  elVec *= area / 3.0;
  return elVec;
#else
  return Eigen::VectorXd::Zero(6);
#endif
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd solveGradientProjection(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    const Eigen::VectorXd& u_h, const lf::assemble::DofHandler& vec_dofh) {
#if SOLUTION
  const auto N_vec_dofs = vec_dofh.NumDofs();

  // Assemble Galerkin mass matrix
  constexpr int cell_codim = 0;
  lf::assemble::COOMatrix<double> M_COO(N_vec_dofs, N_vec_dofs);
  GradientProjectionMassMatrixProvider elMat_builder;
  lf::assemble::AssembleMatrixLocally(cell_codim, vec_dofh, vec_dofh,
                                      elMat_builder, M_COO);

  // Assemble RHS with load grad(u_h)
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_vec_dofs);
  GradientProjectionRhsVectorProvider elVec(fe_space_p, u_h);
  lf::assemble::AssembleVectorLocally(cell_codim, vec_dofh, elVec, phi);

  // Solve M * g = phi
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  return solver.solve(phi);
#else
  return Eigen::VectorXd::Zero(1);
#endif
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::MatrixXd integrateRobotPath(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    const lf::assemble::DofHandler& dofh, const Eigen::VectorXd& sol,
    const lf::assemble::DofHandler& vec_dofh, const Eigen::VectorXd& grad,
    const Eigen::Vector2d& x_start, double dt, int max_steps,
    double u_door_threshold) {
#if SOLUTION
  // Helper: compute barycentric coordinates of point p in triangle
  // with vertices v0, v1, v2. Returns (lambda0, lambda1, lambda2).
  auto barycentric = [](const Eigen::Vector2d& p, const Eigen::Vector2d& v0,
                        const Eigen::Vector2d& v1,
                        const Eigen::Vector2d& v2) -> Eigen::Vector3d {
    Eigen::Matrix2d T;
    T.col(0) = v0 - v2;
    T.col(1) = v1 - v2;
    Eigen::Vector2d lam01 = T.inverse() * (p - v2);
    return {lam01(0), lam01(1), 1.0 - lam01(0) - lam01(1)};
  };

  // Helper: find the triangle containing point p and evaluate
  // the interpolated gradient and potential at p.
  // Returns {gradb_x, gradb_y, u} or {0,0,-1} if not found.
  constexpr int cell_codim = 0;
  constexpr int node_subcodim = 2;
  auto evaluate = [&](const Eigen::Vector2d& p) -> Eigen::Vector3d {
    for (const lf::mesh::Entity* cell : mesh_p->Entities(cell_codim)) {
      if (cell->RefEl() != lf::base::RefEl::kTria()) continue;
      auto corners = lf::geometry::Corners(*(cell->Geometry()));
      Eigen::Vector2d v0 = corners.col(0);
      Eigen::Vector2d v1 = corners.col(1);
      Eigen::Vector2d v2 = corners.col(2);
      Eigen::Vector3d lam = barycentric(p, v0, v1, v2);
      // Check if point is inside triangle (with small tolerance)
      if (lam(0) >= -1e-10 && lam(1) >= -1e-10 && lam(2) >= -1e-10) {
        // Interpolate gradient and potential
        auto scal_idx = dofh.GlobalDofIndices(*cell);
        double u_val = (lam(0) * sol(scal_idx[0])) +
                       (lam(1) * sol(scal_idx[1])) +
                       (lam(2) * sol(scal_idx[2]));
        auto nodes = cell->SubEntities(node_subcodim);
        Eigen::Vector2d grad_val = Eigen::Vector2d::Zero();
        for (int k = 0; k < 3; ++k) {
          auto vi = vec_dofh.GlobalDofIndices(*nodes[k]);
          grad_val(0) += lam(k) * grad(vi[0]);
          grad_val(1) += lam(k) * grad(vi[1]);
        }
        return {grad_val(0), grad_val(1), u_val};
      }
    }
    return {0.0, 0.0, -1.0};  // not found
  };

  // RK4 integration
  std::vector<Eigen::Vector2d> path;
  path.push_back(x_start);
  Eigen::Vector2d x = x_start;

  for (int step = 0; step < max_steps; ++step) {
    // Evaluate -grad(u) at current position
    Eigen::Vector3d ev = evaluate(x);
    if (ev(2) < -0.5) break;              // point outside mesh
    if (ev(2) < u_door_threshold) break;  // reached door

    // RK4 stages: dx/dt = -grad(u)
    Eigen::Vector2d k1 = {-ev(0), -ev(1)};

    ev = evaluate(x + 0.5 * dt * k1);
    if (ev(2) < -0.5) break;
    Eigen::Vector2d k2 = {-ev(0), -ev(1)};

    ev = evaluate(x + 0.5 * dt * k2);
    if (ev(2) < -0.5) break;
    Eigen::Vector2d k3 = {-ev(0), -ev(1)};

    ev = evaluate(x + dt * k3);
    if (ev(2) < -0.5) break;
    Eigen::Vector2d k4 = {-ev(0), -ev(1)};

    x += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    path.push_back(x);
  }

  // Convert to Nx2 matrix
  Eigen::MatrixXd result(path.size(), 2);
  for (size_t i = 0; i < path.size(); ++i) {
    result.row(static_cast<Eigen::Index>(i)) = path[i].transpose();
  }
  return result;
#else
  return Eigen::MatrixXd::Zero(1, 2);
#endif
}
/* SAM_LISTING_END_6 */

}  // namespace RobotNavigation
