/**
 * @file advectionsupg.h
 * @brief NPDE homework AdvectionSUPG code
 * @author R. Hiptmair
 * @date July 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef ADVSUPG_H_H
#define ADVSUPG_H_H

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <algorithm>

namespace AdvectionSUPG {
/** @brief Element matrix provider for streamline-upwind finite-element method
   for for pure advection using quadratic Lagrangian finite elements.
 */
template <class MESHFUNCTION_V>
class SUAdvectionElemMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  SUAdvectionElemMatrixProvider(const SUAdvectionElemMatrixProvider&) = delete;
  SUAdvectionElemMatrixProvider(SUAdvectionElemMatrixProvider&&) noexcept =
      default;
  SUAdvectionElemMatrixProvider& operator=(
      const SUAdvectionElemMatrixProvider&) = delete;
  SUAdvectionElemMatrixProvider& operator=(SUAdvectionElemMatrixProvider&&) =
      delete;
  SUAdvectionElemMatrixProvider(MESHFUNCTION_V& v, bool use_delta = true);
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  ElemMat Eval(const lf::mesh::Entity& cell);
  virtual ~SUAdvectionElemMatrixProvider() = default;

 private:
  // six-point quadrature rule
  const lf::quad::QuadRule qr_{lf::quad::make_TriaQR_P6O4()};
  // A mesh function object providing the velocity field
  MESHFUNCTION_V& v_;
  // Values of reference shape functions at quadrature points
  Eigen::MatrixXd val_ref_lsf_;
  // Gradients of reference shape functions at quadrature points
  Eigen::MatrixXd grad_ref_lsf_;
  // Flag for controlling use of delta scaling
  bool use_delta_;
};

// Constructor
/* SAM_LISTING_BEGIN_1 */
template <class MESHFUNCTION_V>
SUAdvectionElemMatrixProvider<MESHFUNCTION_V>::SUAdvectionElemMatrixProvider(
    MESHFUNCTION_V& v, bool use_delta)
    : v_(v), use_delta_(use_delta) {
  const lf::uscalfe::FeLagrangeO2Tria<double> ref_fe;
  LF_ASSERT_MSG(ref_fe.RefEl() == lf::base::RefEl::kTria(),
                "Only implemented for triangles");
  LF_ASSERT_MSG(ref_fe.NumRefShapeFunctions() == 6,
                "Quadratic Lagrange FE must have six local shape functions");
  LF_ASSERT_MSG(qr_.RefEl() == lf::base::RefEl::kTria(),
                "Quadrature rule must be for triangles");
  LF_ASSERT_MSG(qr_.NumPoints() == 6, "Six-point quadrature rule required");
  // Obtain values and gradients of reference shape functions at
  // quadrature nodes, see \lref{par:lfppparfe}
  val_ref_lsf_ = ref_fe.EvalReferenceShapeFunctions(qr_.Points());
  grad_ref_lsf_ = ref_fe.GradientsReferenceShapeFunctions(qr_.Points());
}
/* SAM_LISTING_END_1 */

// Evaluation operator
/* SAM_LISTING_BEGIN_2 */
template <class MESHFUNCTION_V>
typename SUAdvectionElemMatrixProvider<MESHFUNCTION_V>::ElemMat
SUAdvectionElemMatrixProvider<MESHFUNCTION_V>::Eval(
    const lf::mesh::Entity& cell) {
  LF_ASSERT_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Only implemented for triangles");
  // Obtain geometry information
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr->DimGlobal() == 2,
                "Only implemented for planar triangles");
  // Gram determinant at quadrature points
  const Eigen::VectorXd dets(geo_ptr->IntegrationElement(qr_.Points()));
  LF_ASSERT_MSG(dets.size() == qr_.NumPoints(),
                "Mismatch " << dets.size() << " <-> " << qr_.NumPoints());
  // Fetch the transformation matrices for the gradients
  const Eigen::MatrixXd JinvT(geo_ptr->JacobianInverseGramian(qr_.Points()));
  LF_ASSERT_MSG(JinvT.cols() == 2 * qr_.NumPoints(),
                "Mismatch " << JinvT.cols() << " <-> " << 2 * qr_.NumPoints());
  // Obtain values of velocity at the 6 quadrature points
  std::vector<Eigen::Vector2d> v_vals = v_(cell, qr_.Points());
  // 6 x 6 element matrix (There are six local shape functions.)
  ElemMat mat(6, 6);
  mat.setZero();
  // Compute the factor $\cob{\delta_K}$
  const Eigen::MatrixXd vertices{lf::geometry::Corners(*geo_ptr)};
  // Size of triangle
  const double hK = std::max({(vertices.col(1) - vertices.col(0)).norm(),
                              (vertices.col(2) - vertices.col(1)).norm(),
                              (vertices.col(0) - vertices.col(2)).norm()});
  // Maximal modulus of velocity in quadrature nodes
  const double vmax =
      std::max_element(v_vals.begin(), v_vals.end(),
                       [](Eigen::Vector2d& a, Eigen::Vector2d& b) -> bool {
                         return (a.norm() < b.norm());
                       })
          ->norm();
  const double delta_K = use_delta_ ? (hK / (2.0 * vmax)) : 1.0;
#if SOLUTION
  // Outer summation loop over quadrature points
  for (unsigned int l = 0; l < 6; ++l) {
    // Metric factor $\cob{|\det\Derv\Phibf_K(\wh{\zetabf}_{\ell})}$ scaled with
    // quadrature weight $\cob{\omega_{\ell}}$
    const double fac = qr_.Weights()[l] * dets[l];
    // Compute transformed gradients $\cob{\Vt_{\ell}^i}$ of all local shape
    // functions and collect them in the columns of a 2 x 6-matrix
    const Eigen::Matrix<double, 2, 6> TrfG =
        JinvT.block(0, 2 * l, 2, 2) *
        (grad_ref_lsf_.block(0, 2 * l, 6, 2).transpose());
    // Compute the inner products of the transformed gradients with the
    // velocity vector at the current quadrature point.
    const Eigen::Matrix<double, 1, 6> mvec = v_vals[l].transpose() * TrfG;
    // Loop over columns of element matrix
    for (unsigned int j = 0; j < 6; ++j) {
      // Loop over rows of element matrix
      for (unsigned int i = 0; i < 6; ++i) {
        mat(i, j) +=
            fac * (delta_K * mvec[i] * mvec[j] + mvec[j] * val_ref_lsf_(i, l));
      }
    }
  }
#else
// ========================================
// Your code here
// ========================================
#endif
  return mat;
}
/* SAM_LISTING_END_2 */

/**
 * @brief implements the solution of a pure transport problem , relies on
 * SUAdvectionElemMatrixProvider
 * @param fe_space shared pointer to the finite element space the pde is solved
 * on
 * @return the vector of expansion coefficients of the finite-element solution
 * with respect to the standard nodal basis of S_0^2(M)
 */
Eigen::VectorXd solveRotSUPG(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space);

void cvgL2SUPG();

void visSolution(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space,
    Eigen::VectorXd& u);
}  // namespace AdvectionSUPG

#endif  // #ifndef ADVSUPG_H_H
