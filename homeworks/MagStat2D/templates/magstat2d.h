/**
 * @file magstat2d.h
 * @brief NPDE homework MagStat2D code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef MagStat2D_H_
#define MagStat2D_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembly_types.h>
#include <lf/assemble/fix_dof.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/special_entity_sets.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cstddef>
#include <iostream>
#include <memory>

namespace MagStat2D {
/**
 * @brief Element matrix provider for monolithic Whitney finite element
 * discretization of the 2D magnetostatic variational saddle point problem
 *
 * This class fits the concept of ELEMENT_MATRIX_PROVIDER
 */

/* SAM_LISTING_BEGIN_1 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
class MagStat2DElementMatrixProvider {
 public:
  // The size of the element matrix is $6\times 6$.
  using ElemMat = Eigen::Matrix<double, 6, 6>;
  MagStat2DElementMatrixProvider(const MagStat2DElementMatrixProvider&) =
      delete;
  MagStat2DElementMatrixProvider(MagStat2DElementMatrixProvider&&) noexcept =
      default;
  MagStat2DElementMatrixProvider& operator=(
      const MagStat2DElementMatrixProvider&) = delete;
  MagStat2DElementMatrixProvider& operator=(MagStat2DElementMatrixProvider&&) =
      delete;
  MagStat2DElementMatrixProvider() = delete;
  virtual ~MagStat2DElementMatrixProvider() = default;

  MagStat2DElementMatrixProvider(MESH_FUNCTION mu) : mu_(mu) {}

  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);

 private:
  ElemMat MK_;
  MESH_FUNCTION mu_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
typename MagStat2DElementMatrixProvider<MESH_FUNCTION>::ElemMat
MagStat2DElementMatrixProvider<MESH_FUNCTION>::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  LF_VERIFY_MSG(cell.Geometry()->isAffine(),
                "Triangle must have straight edges");
  // First retrieve value of coefficient $\mu$ at the center of gravity
  // of the triangle
  const Eigen::MatrixXd mpc{
      (Eigen::MatrixXd(2, 1) << 1.0 / 3.0, 1.0 / 3.0).finished()};
  const auto muvals{mu_(cell, mpc)};
  const double muval = muvals[0];
  LF_ASSERT_MSG(muvals.size() == 1, "Only 1 value for mu requested!");
  LF_ASSERT_MSG(muval > 0, "Coefficient must be positive");
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Compute gradients of barycentric coordinate functions and store them in the
  // columns of the $2\times 3$-matrix G
  // clang-format off
  const Eigen::MatrixXd dpt = (Eigen::MatrixXd(2, 1) << 0.0, 0.0).finished();
  const Eigen::Matrix<double, 2, 3> G =
      cell.Geometry()->JacobianInverseGramian(dpt).block(0, 0, 2, 2) *
    (Eigen::Matrix<double, 2, 3>(2,3) << -1, 1, 0,
                                         -1, 0, 1).finished();
  // clang-format on
  // Compute the element matrix  for $-\Delta$ and $\cob{\Cs^0_1}$.
  // See also \lref{mc:ElementMatrixLaplLFE}.
  Eigen::Matrix<double, 3, 3> L = area * G.transpose() * G;
  return MK_;
}
/* SAM_LISTING_END_2 */

/**
 * @brief ENTITY_VECTOR_PROVIDER class for source terms
 * in Whitney FEM for 2D magnetostatic variational saddle point problem
 *
 */
/* SAM_LISTING_BEGIN_7 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
class MagStat2DElementVectorProvider {
 public:
  using ElemVec = Eigen::Matrix<double, 6, 1>;
  MagStat2DElementVectorProvider(const MagStat2DElementVectorProvider&) =
      delete;
  MagStat2DElementVectorProvider(MagStat2DElementVectorProvider&&) noexcept =
      default;
  MagStat2DElementVectorProvider& operator=(
      const MagStat2DElementVectorProvider&) = delete;
  MagStat2DElementVectorProvider& operator=(MagStat2DElementVectorProvider&&) =
      delete;
  virtual ~MagStat2DElementVectorProvider() = default;

  MagStat2DElementVectorProvider(MESH_FUNCTION j_source)
      : j_source_(j_source) {}
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemVec Eval(const lf::mesh::Entity& cell);

 private:
  ElemVec phiK_;
  MESH_FUNCTION j_source_;
};
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
typename MagStat2DElementVectorProvider<MESH_FUNCTION>::ElemVec
MagStat2DElementVectorProvider<MESH_FUNCTION>::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Compute gradients of barycentric coordinate functions, see
  // \lref{cpp:gradbarycoords}. Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns.
  // Note that $\grad\lambda_i$ is accessed as G.col(i-1).
  const auto G{X.inverse().block<2, 3>(1, 0)};
  // Matrix collecting the \textbf{reference coordinates} of the midpoints of
  // the edges of the triangle
  // clang-format off
    const Eigen::MatrixXd mp{(Eigen::MatrixXd(2, 3) <<
			      0.5, 0.5, 0.0,
			      0.0, 0.5, 0.5).finished()};
  // clang-format on
  // Lambda functions in terms of reference coordinates for local Whitney
  // 1-forms, see \prbcref{eq:lsf1}.
  const std::array<std::function<Eigen::Vector2d(Eigen::Vector2d)>, 3> beta{
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return (1.0 - c[0] - c[1]) * G.col(1) - c[0] * G.col(0);
      },
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return c[0] * G.col(2) - c[1] * G.col(1);
      },
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return c[1] * G.col(0) - (1.0 - c[0] - c[1]) * G.col(2);
      }};
  // Obtain values of $\Vf$ at midpoints of edges
  const std::vector<Eigen::Vector2d> f_vals{j_source_(cell, mp)};
  LF_VERIFY_MSG(f_vals.size() == 3, "Too few f values");
  // Fill entries of element vetor. Note that the first three components remain
  // set to zero, because the r.h.s. functional vanishes for the 0-form
  // component of the test function
  phiK_.setZero();
  for (int i = 0; i < 3; ++i) {
    // Iterate over quadrature nodes (midpoints of edges)
    for (int k = 0; k < 3; ++k) {
      phiK_[i + 3] += beta[i](mp.col(k)).dot(f_vals[k]);
    }
  }
  // Multiply with quadrature weight
  phiK_ *= (area / 3.0);
  // Correct for orientation mismatch
  auto relor = cell.RelativeOrientations();
  LF_ASSERT_MSG(relor.size() == 3, "Triangle should have 3 edges!?");
  for (int k = 0; k < 3; ++k) {
    if (relor[k] == lf::mesh::Orientation::negative) {
      // Flip sign of coefficient
      phiK_[k + 3] *= -1;
    }
  }
  return phiK_;
}
/* SAM_LISTING_END_8 */

/**
 * @brief Compute right-hand side vector for 2D magnetostatic
 * Whitney FEM
 *
 * @param dofh DofHandler for monolithic Whitney finite element space for
 * magnetostatic variational saddle point problem
 * @param f MeshFunction providing source vector field
 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
Eigen::VectorXd computeMagStat2DRhsVector(const lf::assemble::DofHandler& dofh,
                                          MESH_FUNCTION j_source) {
  // Total number of FE d.o.f.s
  lf::assemble::size_type N = dofh.NumDofs();
  // Right-hand side vector
  Eigen::VectorXd rhs(N);
  // Object in charge of computing the element vectors
  MagStat2DElementVectorProvider<MESH_FUNCTION> hlevp(j_source);
  // Carry out assembly
  rhs.setZero();
  lf::assemble::AssembleVectorLocally(0, dofh, hlevp, rhs);
  return rhs;
}

/**
 * @brief Assembly of full Galerkin matrix in triplet format
 *
 * @param dofh DofHandler object  for all FE spaces
 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
lf::assemble::COOMatrix<double> buildMagStat2DGalerkinMatrix(
    const lf::assemble::DofHandler& dofh, MESH_FUNCTION mu) {
  // Total number of FE d.o.f.s
  lf::assemble::size_type N = dofh.NumDofs();
  // Full Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A(N, N);
  // Set up computation of element matrix
  MagStat2DElementMatrixProvider hlemp(mu);
  // Assemble \cor{full} Galerkin matrix for Whitney FEM for
  // Hodge-Laplacian mixed variational problem
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, hlemp, A);
  return A;
}

/**
 * @brief Compute Whitney FEM solution of 2D magnetostatic BVP
 *
 * @param dofh DofHandler for monolithic Whitney finite element space for
 * Hodge-Laplacian mixed variational problem
 * @param f MeshFunction providing source vector field
 */
/* SAM_LISTING_BEGIN_3 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION_J,
          lf::mesh::utils::MeshFunction MESH_FUNCTION_MU>
Eigen::VectorXd solveMagStat2DBVP(const lf::assemble::DofHandler& dofh,
                                  MESH_FUNCTION_MU mu,
                                  MESH_FUNCTION_J j_source) {
  // Size of linear system
  const size_t N_dofs = dofh.NumDofs();
  // Galerkin matrix in COO format
  lf::assemble::COOMatrix<double> M_COO{
      buildMagStat2DGalerkinMatrix<MESH_FUNCTION_MU>(dofh, mu)};
  // Right-hand side vector
  Eigen::VectorXd phi{
      computeMagStat2DRhsVector<MESH_FUNCTION_J>(dofh, j_source)};
  // Enforce zero essential boundary conditions
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  // Set boundary flags for both mesh nodes and edges
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
  // Flag all d.o.f.s associated with mesh entities on the boundary
  std::vector<std::pair<bool, double>> ess_dof_select{};
  /* **********************************************************************
     Your code here
     ********************************************************************** */
  // Modify linear system of equations enforcing zero value for d.o.f.s
  // belonging to mesh entities on the boundary
  lf::assemble::FixFlaggedSolutionCompAlt(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      M_COO, phi);

  // Solve linear system using Eigen's sparse direct elimination
  // Examine return status of solver in case the matrix is singular
  const Eigen::SparseMatrix<double> M_crs = M_COO.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  const Eigen::VectorXd dofvec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  return dofvec;
}
/* SAM_LISTING_END_3 */

/** @brief Mesh function providing a proxy vector field for Whitney 1-forms in
 * 2D.
 *
 * Note that the monolithic Whitney FEM d.o.f. layout for the HL mixed
 * variational problem is used. The coefficient vector comprises two parts, the
 * first "node-associated" gives the basis expansion coefficients of the 0-form
 * component, the second the Whitney 1-form expansion coefficients of the 1-form
 * component.
 *
 * Only the second part of the coefficient vector is used by objects of this
 * MeshFunction type.
 *
 */
class MeshFunctionWF1 {
 public:
  MeshFunctionWF1(const lf::assemble::DofHandler& dofh, Eigen::VectorXd coeffs)
      : dofh_(dofh), coeffs_(std::move(coeffs)) {
    const lf::mesh::Mesh& mesh = *dofh.Mesh();
    LF_ASSERT_MSG(dofh_.NumDofs() == coeffs_.size(),
                  "Size mismatch for coeff vector");
    LF_ASSERT_MSG(dofh.NumDofs() == (mesh.NumEntities(2) + mesh.NumEntities(1)),
                  "DofH must manage 1 dof/node and 1 dof/edge");
  }

  // Evaluation operator: returns the values of the vectorfield in the space of
  // Whitney 1-forms at a number of points inside a cell
  std::vector<Eigen::Vector2d> operator()(const lf::mesh::Entity& cell,
                                          const Eigen::MatrixXd& local) const;

 private:
  const lf::assemble::DofHandler& dofh_;
  Eigen::VectorXd coeffs_;
};

/** @brief Mesh function providing a proxy vector field for Whitney 0-forms in
 * 2D.
 *
 * Note that the monolithic Whitney FEM d.o.f. layout for the HL mixed
 * variational problem is used. The coefficient vector comprises two parts, the
 * first "node-associated" gives the basis expansion coefficients of the 0-form
 * component, the second the Whitney 1-form expansion coefficients of the 1-form
 * component.
 *
 * Only the first part of the coefficient vector is used by objects of this
 * MeshFunction type.
 *
 */
class MeshFunctionWF0 {
 public:
  MeshFunctionWF0(const lf::assemble::DofHandler& dofh, Eigen::VectorXd coeffs)
      : dofh_(dofh), coeffs_(std::move(coeffs)) {
    const lf::mesh::Mesh& mesh = *dofh.Mesh();
    LF_ASSERT_MSG(dofh_.NumDofs() == coeffs_.size(),
                  "Size mismatch for coeff vector");
    LF_ASSERT_MSG(dofh.NumDofs() == (mesh.NumEntities(2) + mesh.NumEntities(1)),
                  "DofH must manage 1 dof/node and 1 dof/edge");
  }

  // Evaluation operator: returns the values of a function in the space of
  // Whitney 0-forms at a number of points inside a cell
  std::vector<double> operator()(const lf::mesh::Entity& cell,
                                 const Eigen::MatrixXd& local) const;

 private:
  const lf::assemble::DofHandler& dofh_;
  Eigen::VectorXd coeffs_;
};

/** @brief Approximate "nodal interpolation" of a vectorfield into the
 * Whitney finite element space of 1-forms using Option I 2D Euclidean
 * vector proxies.
 *
 * Special implementation for monolithic handling of d.o.f.s for
 * Hodge-Laplacian mixed variational problem.
 *
 * @param DofHandler for Whitney finite element space of 0-forms and 1-forms
 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
Eigen::VectorXd nodalProjectionWF1(const lf::assemble::DofHandler& dofh,
                                   MESH_FUNCTION vf) {
  const lf::mesh::Mesh& mesh = *dofh.Mesh();
  const Eigen::Index N = dofh.NumDofs();
  LF_ASSERT_MSG(N == (mesh.NumEntities(2) + mesh.NumEntities(1)),
                "DofH must manage 1 dof/node and 1 dof/edge");
  Eigen::VectorXd np_coeffs(N);
  np_coeffs.setZero();
  // Reference coordinates of midpoint of edge
  Eigen::MatrixXd mpc(1, 1);
  mpc(0, 0) = 0.5;
  for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
    // Obtain number of global d.o.f. assciated with edge
    std::span<const lf::assemble::gdof_idx_t> edofs{
        dofh.InteriorGlobalDofIndices(*edge)};
    LF_ASSERT_MSG(edofs.size() == 1, "Only one d.o.f. per edge is allowed");
    const Eigen::Matrix2d endpt = lf::geometry::Corners(*edge->Geometry());
    const Eigen::Vector2d edvec = endpt.col(1) - endpt.col(0);
    // Midpoint quadrature rule:
    const auto vf_vals{vf(*edge, mpc)};
    np_coeffs[edofs[0]] = edvec.dot(vf_vals[0]);
  }
  return np_coeffs;
}

/** @brief test of convergence based on manufactured solution
 *
 * Computes L2 norm of error of the FEM solution for the 1-form component on a
 * sequence of meshes generated by regular refinement.
 */
void testCvgMagStat2DWhitneyFEM(unsigned int refsteps);

}  // namespace MagStat2D

#endif
