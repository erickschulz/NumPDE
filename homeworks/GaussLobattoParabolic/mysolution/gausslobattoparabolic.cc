/**
 * @file gausslobattoparabolic.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "gausslobattoparabolic.h"

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <functional>
#include <memory>
#include <utility>

namespace GaussLobattoParabolic {

/* SAM_LISTING_BEGIN_1 */
lf::assemble::COOMatrix<double> initMbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space) {
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  //====================
  // Your code goes here
  // Replace this dummy assignment for M:
  int N = dofh.NumDofs();
  lf::assemble::COOMatrix<double> M(N, N);
  //====================

  return M;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
lf::assemble::COOMatrix<double> initAbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space) {
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  //====================
  // Your code goes here
  // Replace this dummy assignment for A:
  int N = dofh.NumDofs();
  lf::assemble::COOMatrix<double> A(N, N);
  //====================

  return A;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
RHSProvider::RHSProvider(const lf::assemble::DofHandler& dofh,
                         std::function<double(double)> g)
    : g_(std::move(g)) {
  //====================
  // Your code goes here
  //====================
}

Eigen::VectorXd RHSProvider::operator()(double t) const {
  //====================
  // Your code goes here
  // Replace this dummy return value:
  return Eigen::VectorXd(0);
  //====================
}
/* SAM_LISTING_END_3 */

}  // namespace GaussLobattoParabolic
