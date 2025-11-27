/**
 * @ file norms.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch, Liaowang Huang (refactoring)
 * @ date 22/03/2019, 06/01/2020 (refactoring)
 * @ copyright Developed at ETH Zurich
 */

#include "pointevaluationrhs_norms.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>

namespace PointEvaluationRhs {

/* SAM_LISTING_BEGIN_1 */
double computeL2normLinearFE(const lf::assemble::DofHandler& dofh,
                             const Eigen::VectorXd& mu) {
  double result = 0.0;
#if SOLUTION
  int N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> mass_matrix(N_dofs, N_dofs);
  MassLocalMatrixAssembler my_mat_provider{};
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider,
                                      mass_matrix);
  const Eigen::SparseMatrix<double> mass_mat = mass_matrix.makeSparse();
  result = std::sqrt(mu.dot(mass_mat * mu));
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double computeH1seminormLinearFE(const lf::assemble::DofHandler& dofh,
                                 const Eigen::VectorXd& mu) {
  // calculate stiffness matrix by using the already existing local assembler
  // LinearFELaplaceElementMatrix
  double result = 0.0;
#if SOLUTION
  int N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> stiffness_matrix(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider,
                                      stiffness_matrix);
  const Eigen::SparseMatrix<double> stiffness_mat =
      stiffness_matrix.makeSparse();
  result = std::sqrt(mu.dot(stiffness_mat * mu));
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_2 */

Eigen::MatrixXd MassLocalMatrixAssembler::Eval(const lf::mesh::Entity& entity) {
  Eigen::MatrixXd result;
#if SOLUTION
  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  double volume = lf::geometry::Volume(*geo_ptr);

  if (lf::base::RefEl::kTria() == entity.RefEl()) {
    result.resize(3, 3);
    result.setZero();
    // See Lemma 2.7.5.5 for the derivation of these entries
    double diag = volume / 6.0;
    double non_diag = volume / 12.0;
    result << diag, non_diag, non_diag, non_diag, diag, non_diag, non_diag,
        non_diag, diag;
  } else if (lf::base::RefEl::kQuad() == entity.RefEl()) {
    result.resize(4, 4);
    result.setZero();

    // use quad rule to evaluate all possible combinations of products of
    // basis functions
    lf::uscalfe::FeLagrangeO1Quad<double> quad_element{};
    lf::quad::QuadRule my_quad_rule = lf::quad::make_QuadQR_P4O4();

    // Evaluate the basis functions on the quadrature points
    Eigen::MatrixXd point_eval =
        quad_element.EvalReferenceShapeFunctions(my_quad_rule.Points());

    // Evaluate the integration element at the quadrature points
    Eigen::MatrixXd int_el_eval =
        geo_ptr->IntegrationElement(my_quad_rule.Points());

    // weigh each point according to its weight in the quadrature formula and
    // according to its integration element
    Eigen::MatrixXd point_eval_weighted = point_eval *
                                          my_quad_rule.Weights().asDiagonal() *
                                          int_el_eval.asDiagonal();

    // This product gives us the result of the quadrature rule for all
    // possible combinations of basis functions
    result = point_eval.transpose() * point_eval_weighted;
  } else {
    LF_ASSERT_MSG(
        false, "Function only defined for triangular or quadrilateral cells");
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}

}  // namespace PointEvaluationRhs
