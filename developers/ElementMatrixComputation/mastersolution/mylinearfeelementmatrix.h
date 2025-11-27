/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef MYLINEARFEELEMENTMATRIX_H_
#define MYLINEARFEELEMENTMATRIX_H_

#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace ElementMatrixComputation {

class MyLinearFEElementMatrix {
 public:
  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /*
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a 4x4 matrix, containing the element matrix. The bottom row/column
   *         is not used in the case of a triangle.
   */
  Eigen::Matrix<double, 4, 4> Eval(const lf::mesh::Entity& cell);
};

}  // namespace ElementMatrixComputation

#endif  // MYLINEARFEELEMENTMATRIX_H_
