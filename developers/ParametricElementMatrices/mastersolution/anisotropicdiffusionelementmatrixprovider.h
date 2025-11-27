
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

class AnisotropicDiffusionElementMatrixProvider {
 public:
  /** @brief Constructor storing the vector field of modelling anisotropy */
  AnisotropicDiffusionElementMatrixProvider(
      std::function<Eigen::Vector2d(Eigen::Vector2d)> anisotropy_vec_field)
      : anisotropy_vec_field_(anisotropy_vec_field) {};
  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity&) { return true; }
  /** @brief Main method for computing the element matrix
   * @param cell refers to current cell (triangle or quadrilateral) for which
   * the element matrix is desired. The implementation uses local edge-midpoint
   * quadrature rule. */
  Eigen::MatrixXd Eval(const lf::mesh::Entity& cell);

 private:
  // This vector-valued function of the form d:coords -> vector is used to
  // describe anisotropic material properties. It enters for the diffusion
  // tensor k:coords -> matrix as k(x) = I + d(x)d(x)^T for example in
  // heat conduction models.
  std::function<Eigen::Vector2d(Eigen::Vector2d)> anisotropy_vec_field_;
};  // class AnisotropicDiffusionElementMatrixProvider

}  // namespace ParametricElementMatrices
