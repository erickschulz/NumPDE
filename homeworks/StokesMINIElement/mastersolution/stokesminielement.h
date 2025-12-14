/**
 * @file stokesminielement.h
 * @brief NPDE homework StokesMINIElement code
 * @author Ralf Hiptmair
 * @date June 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef StokesMINIElement_H_
#define StokesMINIElement_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>

namespace StokesMINIElement {  // namespace StokesMINIElement
/**
 * @brief Element matrix provider for simple Stokes FEM
 */
/* SAM_LISTING_BEGIN_1 */
class SimpleFEMElementMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, 9, 9>;
  SimpleFEMElementMatrixProvider(const SimpleFEMElementMatrixProvider&) =
      delete;
  SimpleFEMElementMatrixProvider(SimpleFEMElementMatrixProvider&&) noexcept =
      default;
  SimpleFEMElementMatrixProvider& operator=(
      const SimpleFEMElementMatrixProvider&) = delete;
  SimpleFEMElementMatrixProvider& operator=(SimpleFEMElementMatrixProvider&&) =
      delete;
  SimpleFEMElementMatrixProvider() = default;
  virtual ~SimpleFEMElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);
};
/* SAM_LISTING_END_1 */

/**
 * @brief Element matrix provider for MINI element Stokes FEM
 */
/* SAM_LISTING_BEGIN_5 */
class MINIElementMatrixProvider {
 public:
  using ElemMat = Eigen::Matrix<double, 11, 11>;
  MINIElementMatrixProvider(const MINIElementMatrixProvider&) = delete;
  MINIElementMatrixProvider(MINIElementMatrixProvider&&) noexcept = default;
  MINIElementMatrixProvider& operator=(const MINIElementMatrixProvider&) =
      delete;
  MINIElementMatrixProvider& operator=(MINIElementMatrixProvider&&) = delete;
  MINIElementMatrixProvider() = default;
  virtual ~MINIElementMatrixProvider() = default;
  [[nodiscard]] bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);
};
/* SAM_LISTING_END_5 */

/**
 * @brief Convergence test for simple Stokes FEM
 */
void testCvgSimpleFEM(unsigned int refsteps = 5);

/**
 * @brief Convergence test for MINI Stokes FEM
 */
void testCvgMINIFEM(unsigned int refsteps = 5);

}  // namespace StokesMINIElement

#endif
