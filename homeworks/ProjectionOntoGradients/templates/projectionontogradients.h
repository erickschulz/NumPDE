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

  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================

  // II. Build the (full) right hand side vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce homogeneous dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE using an Eigen solver
  //====================
  // Your code goes here
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_5 */

}  // namespace ProjectionOntoGradients
