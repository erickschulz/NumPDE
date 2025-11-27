/** @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans, Erick Schulz (refactoring)
 * @date 13/03/2019, 19/11/2019 (refactoring)
 * @copyright Developed at ETH Zurich */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <iostream>

namespace ParametricElementMatrices {

class FESourceElemVecProvider {
 public:
  /** @brief Constructor storing the basis expansion vector of the variable
   * coefficient and the finite elements space */
  FESourceElemVecProvider(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      Eigen::VectorXd coeff_expansion)
      : fe_space_(fe_space), coeff_expansion_(coeff_expansion) {}
  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity& cell) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell refers to current cell (triangle or quadrilateral) for which
   * the element veector is desired. The implementation uses local edge-midpoint
   * quadrature rule. */
  Eigen::VectorXd Eval(const lf::mesh::Entity& cell);

 private:
  // Linear first-order lagrangian finite element space
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
  // Finite element basis expansion vector of the coefficient function
  Eigen::VectorXd coeff_expansion_;
};  // class FESourceElemVecProvider

}  // namespace ParametricElementMatrices
