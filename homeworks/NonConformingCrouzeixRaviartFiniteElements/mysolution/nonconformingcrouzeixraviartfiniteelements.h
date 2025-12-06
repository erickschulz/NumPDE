/** @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   16.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich */

#ifndef CRREFERENCEFINITEELEMENT_H
#define CRREFERENCEFINITEELEMENT_H

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/uscalfe/uscalfe.h>

namespace NonConformingCrouzeixRaviartFiniteElements {

class CRReferenceFiniteElement final
    : public lf::fe::ScalarReferenceFiniteElement<double> {
 public:
  lf::base::RefEl RefEl() const override;
  unsigned int Degree() const override;
  lf::assemble::size_type NumRefShapeFunctions() const override;
  lf::assemble::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override;
  lf::assemble::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim, lf::base::sub_idx_t subidx) const override;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override;
  Eigen::MatrixXd EvaluationNodes() const override;
  lf::assemble::size_type NumEvaluationNodes() const override;
  Eigen::Matrix<double, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<double, 1, Eigen::Dynamic>& nodvals) const override;
};

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // CRREFERENCEFINITEELEMENT_H
