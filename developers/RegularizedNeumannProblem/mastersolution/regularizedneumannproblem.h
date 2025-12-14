#ifndef REGULARIZEDNEUMANNPROBLEM_H_
#define REGULARIZEDNEUMANNPROBLEM_H_

/**
 * @file regularizedneumannproblem.h
 * @brief NPDE homework RegularizedNeumannProblem code
 * @author Christian Mitsch, Philippe Peter
 * @date March 2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <memory>
#include <utility>

namespace RegularizedNeumannProblem {

template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE_dropDof(
    const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>
        fe_space,
    const FUNCT_F& f, FUNCT_H& h) {
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  const std::size_t N_dofs = dofh.NumDofs();
  // Right-hand-side vector; don't forget to set to zero initially!
  Eigen::VectorXd rhs_vec(N_dofs);
  rhs_vec.setZero();
  // For returning the sparse Galerkin matrix
  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);

  // ASSEMBLE MATRIX
  lf::assemble::COOMatrix<double> A_aux(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A_aux);

  // ASSEMBLE GLOBAL RHS VECOR
  // Volume part $v\mapsto\int_{\Omega}fv\,\mathrm{d}\Bx$ of the rhs linear
  // functional
  lf::uscalfe::ScalarLoadElementVectorProvider<double, FUNCT_F> my_vec_provider(
      fe_space, f);
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider, rhs_vec);

  // Boundary part $v\mapsto\int_{\partial\Omega}hv\,\mathrm{d}S$ of the rhs
  // linear functional. Only edges on the boundary should be included in
  // assembly, which is achieved by passing a suitable selector predicate to the
  // builder object through a lambda function
  auto bd_edges{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 1)};
  lf::uscalfe::ScalarLoadEdgeVectorProvider my_vec_provider_edge(
      fe_space, h,
      [&bd_edges](const lf::mesh::Entity& edge) { return bd_edges(edge); });
  // co-dimension 1 because we locally assemble on edges!
  lf::assemble::AssembleVectorLocally(1, dofh, my_vec_provider_edge, rhs_vec);

  // Now fix the solution to be 0 at the node p
#if SOLUTION
  // Tell the function FixFlaggedSolutionComponents that is should
  // modify the matrix and the right-hand-side vector in a way that amounts to
  // dropping the first tent function from the basis.
  /* SAM_LISTING_BEGIN_1 */
  auto selector = [](lf::base::glb_idx_t idx) -> std::pair<bool, double> {
    if (idx == 0) {
      return {true, 0.0};  // fix first d.o.f. to zero
    }
    return {false, 42.0};  // keep all others
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_aux, rhs_vec);
  /* SAM_LISTING_END_1 */
#else
  //====================
  // Your code goes here
  //====================
#endif
  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

// You can write a helper class which should implement
// ENTITY_VECTOR_PROVIDER to calculate vector c using LehrFEM assembly functions
/* SAM_LISTING_BEGIN_6 */
class VecHelper {
 public:
  explicit VecHelper() {}
  bool isActive(const lf::mesh::Entity& entity) const { return true; }
  Eigen::Vector3d Eval(const lf::mesh::Entity& entity) {
    LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                  "Function only defined for triangular cells");
    Eigen::Vector3d result;
#if SOLUTION
    // Obtain shape information for the cell
    const lf::geometry::Geometry* geo_ptr = entity.Geometry();
    // Fetch area
    const double area = lf::geometry::Volume(*geo_ptr);
    // Initialize element vector |K|/3*[1,1,1]^T
    result = (area / 3.0) * Eigen::Vector3d::Ones();
#else
    //====================
    // Your code goes here
    //====================
#endif
    return result;
  }
};
/* SAM_LISTING_END_6 */

Eigen::VectorXd assembleVector_c(const lf::assemble::DofHandler& dofh);

template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE_augment(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const FUNCT_F& f, const FUNCT_H& h) {
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  const std::size_t N_dofs = dofh.NumDofs() + 1;
  Eigen::VectorXd rhs_vec(N_dofs);
  rhs_vec.setZero();
  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);

  // ASSEMBLE MATRIX
  lf::assemble::COOMatrix<double> A_aux(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A_aux);

  // ASSEMBLE GLOBAL RHS VECOR
  // Volume part $v\mapsto\int_{\Omega}fv\,\mathrm{d}\Bx$ of the rhs linear
  // functional
  lf::uscalfe::ScalarLoadElementVectorProvider<double, FUNCT_F> my_vec_provider(
      fe_space, f);
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider, rhs_vec);

  // Boundary part $v\mapsto\int_{\partial\Omega}hv\,\mathrm{d}S$ of the rhs
  // linear functional. Only edges on the boundary should be included in
  // assembly, which is achieved by passing a suitable selector predicate to the
  // builder object through a lambda function
  auto bd_edges{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 1)};
  lf::uscalfe::ScalarLoadEdgeVectorProvider my_vec_provider_edge(
      fe_space, h,
      [&bd_edges](const lf::mesh::Entity& edge) { return bd_edges(edge); });
  // co-dimension 1 because we locally assemble on edges
  lf::assemble::AssembleVectorLocally(1, dofh, my_vec_provider_edge, rhs_vec);

  // Calculate c and insert it into A
#if SOLUTION
  /* SAM_LISTING_BEGIN_2 */
  // Calculate vector c
  Eigen::VectorXd c = assembleVector_c(dofh);
  // Add c to the matrix
  for (int i = 0; i < c.size(); i++) {
    // add triplets for vector c
    A_aux.AddToEntry(N_dofs - 1, i, c(i));
    A_aux.AddToEntry(i, N_dofs - 1, c(i));
  }
  /* SAM_LISTING_END_2 */
#else
  //====================
  // Your code goes here
  //====================
#endif

  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

}  // namespace RegularizedNeumannProblem

#endif  // REGULARIZEDNEUMANNPROBLEM_H_
