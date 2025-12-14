#ifndef STANDARD_FEM_H
#define STANDARD_FEM_H

/**
 * @file standard_fem.h
 * @brief Solves the CD BVP based on a standard FEM approach
 * @author Philippe Peter
 * @date June 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <memory>

#include "convection_emp.h"

namespace ConvectionDiffusion {

/**
 * @brief Solves the Convection-Diffusion BVP with nonhomogeneous Dirichlet
 * boundary conditions using a standard FEM approach
 */
template <typename DIFFUSION_COEFF, typename CONVECTION_COEFF,
          typename FUNCTOR_F, typename FUNCTOR_G>
Eigen::VectorXd SolveCDBVPStandardFem(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space,
    DIFFUSION_COEFF eps, CONVECTION_COEFF v, FUNCTOR_F f, FUNCTOR_G g) {
  // Wrap functions into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};
  lf::mesh::utils::MeshFunctionGlobal mf_eps{eps};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // mesh and dofhanlder
  auto mesh_p = fe_space->Mesh();
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

  // ASSEMBLE GALERKIN MATRIX
  // First the part corresponding to the laplacian
  lf::uscalfe::ReactionDiffusionElementMatrixProvider laplacian_provider(
      fe_space, mf_eps, lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider, A);

  // Next part corresponding to the convection term:
  ConvectionElementMatrixProvider convection_provider(v);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider, A);

  // RIGHT-HAND SIDE VECTOR
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();
  lf::uscalfe::ScalarLoadElementVectorProvider elvec_provider(fe_space, mf_f);
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_provider, phi);

  // IMPOSE DIRICHLET BC
  // Obtain specification for shape functions on edges
  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Create a dataset of boolean flags indicating edges on the boundary of the
  // mesh
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // boundary.
  auto ess_bdc_flags_values{
      lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g)};

  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::uscalfe::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  // SOLVE LINEAR SYSTEM
  Eigen::SparseMatrix A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  return sol_vec;
}

}  // namespace ConvectionDiffusion

#endif  // STANDARD_GALERKIN_H