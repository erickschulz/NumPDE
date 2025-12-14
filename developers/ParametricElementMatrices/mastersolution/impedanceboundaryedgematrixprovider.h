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
#include <stdexcept>

namespace ParametricElementMatrices {

class ImpedanceBoundaryEdgeMatrixProvider {
 public:
  /** @brief Constructor storing the basis expansion vector of the variable
   * coefficient, the finite elements space and the boundary edge predicate */
  ImpedanceBoundaryEdgeMatrixProvider(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      Eigen::VectorXd coeff_expansion);
  bool isActive(const lf::mesh::Entity& edge);
  Eigen::MatrixXd Eval(const lf::mesh::Entity& edge);

 private:
  // Linear first-order lagrangian finite element space
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
  // Finite element basis expansion vector of the coefficient function
  Eigen::VectorXd coeff_expansion_;
  // Predicate returning true if an edge is on the boundary
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<bool>> bd_flags_;
};
}  // namespace ParametricElementMatrices
