#ifndef BOUNDARYWAVE_HPP
#define BOUNDARYWAVE_HPP

/** @file
 * @brief NPDE BoundaryWave
 * @author Erick Schulz
 * @date 24/07/2019
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/SparseLU>

namespace BoundaryWave {

// Library functions
lf::assemble::COOMatrix<double> buildM(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p);

lf::assemble::COOMatrix<double> buildA(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p);

/* SAM_LISTING_BEGIN_7 */
template <typename FUNCTOR_U, typename FUNCTOR_V>
std::pair<Eigen::VectorXd, Eigen::VectorXd> interpolateInitialData(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNCTOR_U&& u0, FUNCTOR_V&& v0) {
  Eigen::VectorXd dof_vector_u0, dof_vector_v0;

  // Generate Lehrfem++ mesh functions out of the functors
  //====================
  // Your code goes here
  //====================

  std::pair<Eigen::VectorXd, Eigen::VectorXd> initialData =
      std::make_pair(dof_vector_u0, dof_vector_v0);
  return initialData;
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
template <typename FUNCTOR_U, typename FUNCTOR_V>
Eigen::VectorXd solveBoundaryWave(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    FUNCTOR_U&& u0, FUNCTOR_V&& v0, double T, unsigned int N) {
  Eigen::VectorXd bdyWaveSol;

  double step_size = T / N;
  // Obtain initial data
  std::pair<Eigen::VectorXd, Eigen::VectorXd> initialData =
      interpolateInitialData<std::function<double(Eigen::Vector2d)>,
                             std::function<double(Eigen::Vector2d)>>(
          fe_space_p, std::forward<FUNCTOR_U>(u0), std::forward<FUNCTOR_V>(v0));
  // Obtain Galerkin matrices
  lf::assemble::COOMatrix<double> M = buildM(fe_space_p);
  lf::assemble::COOMatrix<double> A = buildA(fe_space_p);
  //====================
  // Your code goes here
  //====================
  return bdyWaveSol;
};
/* SAM_LISTING_END_8 */

}  // namespace BoundaryWave

#endif
