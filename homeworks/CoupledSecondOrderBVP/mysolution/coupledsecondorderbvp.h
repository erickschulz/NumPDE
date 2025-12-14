/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP
 * @author Erick Schulz
 * @date 13/11/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <memory>
#include <stdexcept>
#include <vector>

namespace CoupledSecondOrderBVP {

template <typename SCALAR>
class FeSpaceLagrangeO2 : public lf::uscalfe::UniformScalarFESpace<SCALAR> {
 public:
  FeSpaceLagrangeO2() = delete;
  FeSpaceLagrangeO2(const FeSpaceLagrangeO2&) = delete;
  FeSpaceLagrangeO2(FeSpaceLagrangeO2&&) noexcept = default;
  FeSpaceLagrangeO2& operator=(const FeSpaceLagrangeO2&) = delete;
  FeSpaceLagrangeO2& operator=(FeSpaceLagrangeO2&&) noexcept = default;
  explicit FeSpaceLagrangeO2(
      const std::shared_ptr<const lf::mesh::Mesh>& mesh_p)
      : lf::uscalfe::UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<lf::uscalfe::FeLagrangeO2Tria<SCALAR>>(),
            std::make_shared<lf::uscalfe::FeLagrangeO2Quad<SCALAR>>(),
            std::make_shared<lf::uscalfe::FeLagrangeO2Segment<SCALAR>>(),
            std::make_shared<lf::fe::FePoint<SCALAR>>(2)) {}
  ~FeSpaceLagrangeO2() override = default;
};  // FeSpaceLagrangeO2

/** @Brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin stiffness and mass matrices. It transforms every columns and rows
 * associated to a global index belonging to a degree of freedom lying on the
 * boundary to zero entries but the diagonal one which is set to 1.0
 * @param selectvals is the predicate identifying the boundary indices of the
 * rows and columns that are to be dropped */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsAndColumns(SELECTOR&& selectvals,
                              lf::assemble::COOMatrix<SCALAR>& A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/** @Brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin stiffness and mass matrices. It transforms every rows
 * associated to a global index belonging to a degree of freedom lying on the
 * boundary to zero entries
 * @param selectvals is the predicate identifying the boundary indices of the
 * rows that are to be dropped */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRows(SELECTOR&& selectvals, lf::assemble::COOMatrix<SCALAR>& M) {
  M.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i));
      });
}

// Function solving the coupled BVP
template <typename FUNCTOR>
Eigen::VectorXd solveCoupledBVP(
    std::shared_ptr<FeSpaceLagrangeO2<double>>& fe_space, double gamma,
    FUNCTOR&& f) {
  Eigen::VectorXd sol_vec;  // solution vector
  // Get pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Data related to Dirichlet B.C.
  // Obtain an array of boolean flags for the nodes of the mesh, 'true'
  // indicates that the node lies on the boundary
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Similarly for edges
  auto edges_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsAndColumns
  auto bd_selector = [&nodes_bd_flags, &edges_bd_flags,
                      &dofh](unsigned int idx) -> bool {
    if (dofh.Entity(idx).RefEl() == lf::base::RefElType::kPoint) {
      return nodes_bd_flags(dofh.Entity(idx));
    }
    return edges_bd_flags(dofh.Entity(idx));
  };

  /* I : Creating coefficients as Lehrfem++ mesh functions */
  // Coefficients used in the class template
  // ReactionDiffusionElementMatrixProvider<SCALAR,DIFF_COEFF,REACTION_COEFF>
  auto const_one = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  auto const_zero = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  /* II: Instantiating finite element matrices and right hand side vector*/
  // Matrices in triplet format holding  the Galerkin matrices
  // (set to zero initially)
  lf::assemble::COOMatrix<double> A0(N_dofs, N_dofs);  // upper left block
  lf::assemble::COOMatrix<double> A1(N_dofs, N_dofs);  // lower right block
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);   // off-diag blocks
  lf::assemble::COOMatrix<double> L(2 * N_dofs, 2 * N_dofs);  // full matrix
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);       // source f
  phi.setZero();  // set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhs(2 * N_dofs);  // full rhs vector
  rhs.setZero();  // set to zero initially

  /* III : Computing the Galerkin matrices */
  // III.i Computing A0 : standard negative Laplace Galerkin Matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_zero)>
      A0_builder(fe_space, const_one, const_zero);
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A0_builder, A0);
  // Enforce Dirichlet boundary conditions
  dropMatrixRowsAndColumns(bd_selector, A0);
  // III.ii Computing A01 : Laplacian with reaction term
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_one)>
      A1_builder(fe_space, const_one, const_one);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A1_builder, A1);
  // III.iii Computing M : standard mass matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_zero), decltype(const_one)>
      M_builder(fe_space, const_zero, const_one);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_builder, M);
  // Enforce Dirichlet boundary conditions
  dropMatrixRows(bd_selector, M);
  // IV : Computing the element vector phi (associated to source f)
  // Wrap the lambda source function f in a Lehrfem++ MeshFunction
  auto mf_f = lf::mesh::utils::MeshFunctionGlobal(f);
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);
  // Assigning zero to the boundary values of phi
  for (const lf::mesh::Entity* node : mesh_p->Entities(2)) {
    if (nodes_bd_flags(*node)) {
      auto dof_idx = dofh.GlobalDofIndices(*node);
      LF_ASSERT_MSG(
          dofh.NumLocalDofs(*node) == 1,
          "Too many global indices were returned for a vertex entity!");
      phi(dof_idx[0]) = 0.0;
    }
  }

  /* V : Assemble the full linear system matrix and right hand side vector */
  //                        _        _
  //       L (u  p)^T  :=  |  A0    M | (u  p)^T  = (f 0)^T           (*)
  //                       |_ M^T  A1_|
  //
  // V.i Inserting A0 in L
  const std::vector<Eigen::Triplet<double>> A0_triplets_vec = A0.triplets();
  for (auto& triplet : A0_triplets_vec) {
    L.AddToEntry(triplet.row(), triplet.col(), triplet.value());
  }
  // V.ii Inserting A1 in L
  const std::vector<Eigen::Triplet<double>> A1_triplets_vec = A1.triplets();
  for (auto& triplet : A1_triplets_vec) {
    L.AddToEntry(triplet.row() + N_dofs, triplet.col() + N_dofs,
                 triplet.value());
  }
  // V.iii Inserting M in L
  // Notice the symmetry between the entries. The two indices of the entry
  // are flipped so that it is the tranpose of M that is added to the left lower
  // diagonal block of L.
  const std::vector<Eigen::Triplet<double>> M_triplets_vec = M.triplets();
  for (auto& triplet : M_triplets_vec) {
    L.AddToEntry(triplet.row(), triplet.col() + N_dofs,
                 triplet.value());  // for M in upper right block
    L.AddToEntry(triplet.col() + N_dofs, triplet.row(),
                 triplet.value());  // for M^T in lower left block
  }
  // V.iv Assembling full right hand side vector
  // Recall that rhs was initialized with zero entries
  for (int i = 0; i < N_dofs; i++) {
    rhs(i) = phi(i);
  }

  // VI. Solve the LSE (*) of step (V) using an Eigen solver
  // Convert COO matrix A into CRS format using Eigen's internal
  // conversion routines.
  Eigen::SparseMatrix<double> L_crs = L.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  solver.compute(L_crs);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  sol_vec = solver.solve(rhs);

  return sol_vec;
}

}  // namespace CoupledSecondOrderBVP
