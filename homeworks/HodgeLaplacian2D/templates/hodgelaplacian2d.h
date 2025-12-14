/**
 * @file hodgelaplacian2d.h
 * @brief NPDE homework HodgeLaplacian2D code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef HodgeLaplacian2D_H_
#define HodgeLaplacian2D_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/assembly_types.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace HodgeLaplacian2D {
/**
 * @brief Element matrix provider for monolithic Whitney finite element
 * discretization of the Hodge Laplacian generalized saddle point problem.
 *
 * This class fits the concept of ELEMENT_MATRIX_PROVIDER
 */

/* SAM_LISTING_BEGIN_1 */
class HodgeLaplacian2DElementMatrixProvider {
 public:
  // The size of the element matrix is $6\times 6$.
  using ElemMat = Eigen::Matrix<double, 6, 6>;
  HodgeLaplacian2DElementMatrixProvider(
      const HodgeLaplacian2DElementMatrixProvider&) = delete;
  HodgeLaplacian2DElementMatrixProvider(
      HodgeLaplacian2DElementMatrixProvider&&) noexcept = default;
  HodgeLaplacian2DElementMatrixProvider& operator=(
      const HodgeLaplacian2DElementMatrixProvider&) = delete;
  HodgeLaplacian2DElementMatrixProvider& operator=(
      HodgeLaplacian2DElementMatrixProvider&&) = delete;
  HodgeLaplacian2DElementMatrixProvider() = default;
  virtual ~HodgeLaplacian2DElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  // Computation of element matrix $\VM_K$: two versions
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);
  // Alternative implementation created for debugging purposes
  [[nodiscard]] ElemMat Eval_ref(const lf::mesh::Entity& cell);

 private:
  ElemMat MK_;
};
/* SAM_LISTING_END_1 */

/**
 * @brief ENTITY_VECTOR_PROVIDER class for source terms
 * in Whitney FEM for Hodge-Laplacian mixed variational problem
 *
 */
/* SAM_LISTING_BEGIN_7 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
class HodgeLaplacian2DElementVectorProvider {
 public:
  using ElemVec = Eigen::Matrix<double, 6, 1>;
  HodgeLaplacian2DElementVectorProvider(
      const HodgeLaplacian2DElementVectorProvider&) = delete;
  HodgeLaplacian2DElementVectorProvider(
      HodgeLaplacian2DElementVectorProvider&&) noexcept = default;
  HodgeLaplacian2DElementVectorProvider& operator=(
      const HodgeLaplacian2DElementVectorProvider&) = delete;
  HodgeLaplacian2DElementVectorProvider& operator=(
      HodgeLaplacian2DElementVectorProvider&&) = delete;
  virtual ~HodgeLaplacian2DElementVectorProvider() = default;

  HodgeLaplacian2DElementVectorProvider(MESH_FUNCTION f) : f_(f) {}
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemVec Eval(const lf::mesh::Entity& cell);

 private:
  ElemVec phiK_;
  MESH_FUNCTION f_;
};
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
typename HodgeLaplacian2DElementVectorProvider<MESH_FUNCTION>::ElemVec
HodgeLaplacian2DElementVectorProvider<MESH_FUNCTION>::Eval(
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
  /* **********************************************************************
     Your code here
     ********************************************************************** */
  return phiK_;
}
/* SAM_LISTING_END_8 */

/**
 * @brief Compute right-hand side vector for Hodge-Laplacian
 * Whitney FEM
 *
 * @param dofh DofHandler for monolithic Whitney finite element space for
 * Hodge-Laplacian mixed variational problem
 * @param f MeshFunction providing source vector field
 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
Eigen::VectorXd computeHodgeLaplaceRhsVector(
    const lf::assemble::DofHandler& dofh, MESH_FUNCTION f) {
  // Total number of FE d.o.f.s
  lf::assemble::size_type N = dofh.NumDofs();
  // Right-hand side vector
  Eigen::VectorXd rhs(N);
  // Object in charge of computing the element vectors
  HodgeLaplacian2DElementVectorProvider<MESH_FUNCTION> hlevp(f);
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
lf::assemble::COOMatrix<double> buildHodgeLaplacianGalerkinMatrix(
    const lf::assemble::DofHandler& dofh);

/**
 * @brief Compute Whitney FEM solution of Hodge Laplace BVP
 *
 * @param dofh DofHandler for monolithic Whitney finite element space for
 * Hodge-Laplacian mixed variational problem
 * @param f MeshFunction providing source vector field
 */
template <lf::mesh::utils::MeshFunction MESH_FUNCTION>
Eigen::VectorXd solveHodgeLaplaceBVP(const lf::assemble::DofHandler& dofh,
                                     MESH_FUNCTION f) {
  // Galerkin matrix in COO format
  lf::assemble::COOMatrix<double> M_COO{
      buildHodgeLaplacianGalerkinMatrix(dofh)};
  const Eigen::SparseMatrix<double> M_crs = M_COO.makeSparse();
  // Right-hand side vector
  const Eigen::VectorXd phi{
      computeHodgeLaplaceRhsVector<MESH_FUNCTION>(dofh, f)};

  // Solve linear system using Eigen's sparse direct elimination
  // Examine return status of solver in case the matrix is singular
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  const Eigen::VectorXd dofvec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return dofvec;
}

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
  MeshFunctionWF0(lf::assemble::DofHandler& dofh, Eigen::VectorXd coeffs)
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
  lf::assemble::DofHandler& dofh_;
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

/** @brief Aggregate vector proxy of 1-form component into mesh nodes
 *
 */
std::pair<lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d>>
reconstructNodalFields(const lf::assemble::DofHandler& dofh,
                       const Eigen::VectorXd& coeffs);

/** @brief test of convergence based on manufactured solution
 *
 * Computes L2 norm of error of the FEM solution for the 1-form component on a
 * sequence of meshes generated by regular refinement.
 */
void testCvgHLWhitneyFEM(unsigned int refsteps);

}  // namespace HodgeLaplacian2D

#endif
