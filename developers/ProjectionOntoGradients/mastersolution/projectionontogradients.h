/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author Erick Schulz, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <utility>

namespace ProjectionOntoGradients {

/* SAM_LISTING_BEGIN_1 */
class ElementMatrixProvider {
 public:
  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity& /*entity*/) const { return true; }
};
/* SAM_LISTING_END_1 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
class GradProjRhsProvider {
 public:
  explicit GradProjRhsProvider(FUNCTOR f) : f_(f) {}

  Eigen::Vector3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity& /*entity*/) const { return true; }

 private:
  FUNCTOR f_;
};
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
Eigen::Vector3d GradProjRhsProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  Eigen::Vector3d loc_vec;

#if SOLUTION
  // get area of the entity
  const double area = lf::geometry::Volume(*geo_ptr);

  // calculate center of mass
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const Eigen::Vector2d c = corners.rowwise().sum() / 3.;

  // value of the functor f at center of mass of the entity
  const Eigen::Vector2d func_value = f_(c);

  // calculate the gradients of the basis functions
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // vector of the gradients of the basis function evaluated at c
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  loc_vec = area * (grad_basis.transpose() * func_value);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return loc_vec;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
template <typename FUNCTOR>
Eigen::VectorXd projectOntoGradients(const lf::assemble::DofHandler& dofh,
                                     FUNCTOR f) {
  const lf::assemble::size_type N_dofs = dofh.NumDofs();
  Eigen::VectorXd sol_vec;

  // I. Build the (full) Galerkin matrix for the linear system.
#if SOLUTION
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  ProjectionOntoGradients::ElementMatrixProvider my_mat_provider;
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // II. Build the (full) right hand side vector
#if SOLUTION
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  ProjectionOntoGradients::GradProjRhsProvider my_vec_provider(f);
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider, phi);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // III. Enforce homogeneous dirichlet boundary conditions
#if SOLUTION
  // We do this by selecting the DOFs on the boundary and setting the
  // values to zero. Note, that we could also use this to set the boundary
  // to any other value (if bc were not homogeneous)

  // To select the right DOFs we need a selector:
  // * for all DOFs on the boundary return true and the value on the boundary
  // * for all other DOFs return false (and whatever as value)
  const double boundary_val = 0;
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2)};
  auto my_selector = [&dofh, &bd_flags, &boundary_val](unsigned int dof_idx) {
    if (bd_flags(dofh.Entity(dof_idx))) {
      return std::make_pair(true, boundary_val);
    }  // interior node: the value we return here does not matter
    return std::make_pair(false, 42.0);
  };
  // Since we know the values on the boundary we know the solution on these
  // DOFs and we can write the Galerkin LSE in block format and solve only
  // for the unknown coefficients. This modification is done by the following
  // function FixFlaggedSolutionComponents(). We use the selector we have
  // defined above.
  // See \lref{sec:essbdc} for explanations.
  lf::assemble::FixFlaggedSolutionComponents<double>(my_selector, A, phi);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // IV. Solve the LSE using an Eigen solver
#if SOLUTION
  // Convert from triplet to CRS format
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  sol_vec = solver.solve(phi);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return sol_vec;
}
/* SAM_LISTING_END_5 */

}  // namespace ProjectionOntoGradients
