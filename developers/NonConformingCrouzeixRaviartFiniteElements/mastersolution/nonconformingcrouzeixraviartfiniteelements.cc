/** @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 *  @author Anian Ruoss, edited Amélie Loher
 *  @date   16.03.2019, 03.03.20
 *  @copyright Developed at ETH Zurich */

#include "nonconformingcrouzeixraviartfiniteelements.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

// Crouzeix-Raviart finite element space defined on triangular meshes only
/* SAM_LISTING_BEGIN_1 */
lf::base::RefEl CRReferenceFiniteElement::RefEl() const {
  lf::base::RefElType ref_el_type;
// TODO: task 2-14.q)
#if SOLUTION
  ref_el_type = lf::base::RefElType::kTria;
#else
//====================
// Your code goes here
//====================
#endif
  return lf::base::RefEl(ref_el_type);
}
/* SAM_LISTING_END_1 */

// Crouzeix-Raviart are piecewise linear polynomials
/* SAM_LISTING_BEGIN_2 */
unsigned int CRReferenceFiniteElement::Degree() const {
  unsigned int degree;
// TODO: task 2-14.q)
#if SOLUTION
  degree = 1;
#else
//====================
// Your code goes here
//====================
#endif
  return degree;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions() const {
  lf::assemble::size_type num_ref_shape_functions;
// TODO: task 2-14.q)
#if SOLUTION
  num_ref_shape_functions = 3;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return num_ref_shape_functions;
}

lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions(
    lf::assemble::dim_t codim) const {
  switch (codim) {
    // TODO: task 2-14.q)
#if SOLUTION
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return 0;
#else
//====================
// Your code goes here
//====================
#endif
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}

lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions(
    lf::assemble::dim_t codim, lf::base::sub_idx_t subidx) const {
  switch (codim) {
    // TODO: task 2-14.q)
#if SOLUTION
    case 0:
      LF_VERIFY_MSG((0 == subidx),
                    "Index of cell is out of range for triangle");
      return 0;
    case 1:
      LF_VERIFY_MSG((0 <= subidx && subidx < 3),
                    "Index of edge is out of range for triangle");
      return 1;
    case 2:
      LF_VERIFY_MSG((0 <= subidx && subidx < 3),
                    "Index of vertex is out of range for triangle");
      return 0;
#else
//====================
// Your code goes here
//====================
#endif
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::EvalReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  // Data
  const auto num_points =
      static_cast<lf::assemble::size_type>(refcoords.cols());
  // Initialize a matrix that will store the values of the reference basis
  // functions evaluated at the coordinates passed as arguments
  Eigen::MatrixXd eval_ref_shape_functions(3, num_points);
// TODO: task 2-14.q)
#if SOLUTION
  Eigen::MatrixXd ones = Eigen::VectorXd::Ones(num_points).transpose();
  // Evaluate the basis functions
  eval_ref_shape_functions.row(0) = ones - 2. * refcoords.row(1);
  eval_ref_shape_functions.row(1) = 2. * refcoords.colwise().sum() - ones;
  eval_ref_shape_functions.row(2) = ones - 2. * refcoords.row(0);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return eval_ref_shape_functions;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::GradientsReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  // Data
  const auto num_points =
      static_cast<lf::assemble::size_type>(refcoords.cols());
  // Initialize a matrix that will store the gradients of the reference basis
  // functions evaluated at the coordinates passed as arguments
  Eigen::MatrixXd grad_ref_shape_functions(3, 2 * num_points);
// TODO: task 2-14.r)
#if SOLUTION
  // Evaluate the gradients
  grad_ref_shape_functions.row(0) = (Eigen::Vector2d() << 0, -2)
                                        .finished()
                                        .transpose()
                                        .replicate(1, num_points);
  grad_ref_shape_functions.row(1) =
      2. * Eigen::VectorXd::Ones(2 * num_points).transpose();
  grad_ref_shape_functions.row(2) = (Eigen::Vector2d() << -2, 0)
                                        .finished()
                                        .transpose()
                                        .replicate(1, num_points);
#else
//====================
// Your code goes here
//====================
#endif
  return grad_ref_shape_functions;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::MatrixXd CRReferenceFiniteElement::EvaluationNodes() const {
  Eigen::MatrixXd eval_nodes(2, 3);
// TODO: task 2-14.s)
#if SOLUTION
  eval_nodes << .5, .5, 0, 0, .5, .5;
#else
//====================
// Your code goes here
//====================
#endif
  return eval_nodes;
}
/* SAM_LISTING_END_6 */

lf::assemble::size_type CRReferenceFiniteElement::NumEvaluationNodes() const {
  return 3;
}

/* SAM_LISTING_BEGIN_7 */
Eigen::Matrix<double, 1, Eigen::Dynamic>
CRReferenceFiniteElement::NodalValuesToDofs(
    const Eigen::Matrix<double, 1, Eigen::Dynamic>& nodvals) const {
  LF_VERIFY_MSG(nodvals.cols() == NumEvaluationNodes(),
                "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());

  Eigen::MatrixXd coeffs;
// TODO: task 2-14.s)
#if SOLUTION
  // Linear mapping is identity since the set of reference shape functions
  // forms a cardinal basis with respect to the interpolation nodes
  coeffs = nodvals;
#else
//====================
// Your code goes here
//====================
#endif
  return coeffs;
}
/* SAM_LISTING_END_7 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements
