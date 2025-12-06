/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef MYLINEARLOADVECTOR_H_
#define MYLINEARLOADVECTOR_H_

#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <functional>
#include <utility>

namespace ElementMatrixComputation {

class MyLinearLoadVector {
 public:
  /** @brief Constructor storing the right hand side function */
  explicit MyLinearLoadVector(std::function<double(const Eigen::Vector2d&)> f)
      : f_(std::move(f)) {}

  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /**
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   * The implementation uses simple edge midpoint based quadrature and an
   * approximation of the volume of a cell just using the integration element at
   * the barycenter.
   */
  Eigen::Vector4d Eval(const lf::mesh::Entity& cell);

 private:
  /** `f_(x)` where `x` is a 2D vector that provides the evaluation of the
   * source function */
  std::function<double(const Eigen::Vector2d&)> f_;
};

}  // namespace ElementMatrixComputation

#endif  // MYLINEARLOADVECTOR_H_
