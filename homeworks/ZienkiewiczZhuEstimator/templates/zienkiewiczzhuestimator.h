#ifndef ZIENKIEWICZZHUESTIMATOR_HPP
#define ZIENKIEWICZZHUESTIMATOR_HPP

/** @file
 * @brief NPDE ZienkiewiczZhuEstimator
 * @author Erick Schulz
 * @date 25/07/2019
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ZienkiewiczZhuEstimator {

/* SIMPLIFICATION OF TYPE NAMES */
using coord_t = Eigen::Vector2d;

/* LIBRARY CLASSES */
/** @brief This class implements a Lehrfem++ matrix provider defining a
VectorProjectionMatrixProvider::Eval function returning the local MASS matrix
for linear first-order lagrange FE bases over mixed triangular and quadrilateral
meshes. Integration is performed using appropriate quadrature rule.*/
class VectorProjectionMatrixProvider {
 public:
  /** @brief default constructor */
  explicit VectorProjectionMatrixProvider() = default;
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /** @brief Main method for computing the element matrix
   * @param cell refers to current cell for which the element matrix is desired
   * The implementation uses appropriate quadrature rule of the cells*/
  Eigen::MatrixXd Eval(const lf::mesh::Entity& entity);
};  // class VectorProjectionMatrixProvider

/** @brief This class implements a Lehrfem++ matrix provider defining a
GradientProjectionVectorProvider::Eval function returning the local
contribution to the element vectors for linear first-order lagrange FE bases
over triangular mesh (only!). Integration over the triangular cells is performed
using the trapezoidal rule.*/
class GradientProjectionVectorProvider {
 public:
  /** @brief Constructor storing the right hand side function */
  explicit GradientProjectionVectorProvider(
      const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
      const Eigen::VectorXd& mu)
      : _mu(mu), _fe_space_p(fe_space_p) {}
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell current cell for which the element vector is desired
   * The implementation uses an appropriate quadrature rule.*/
  Eigen::VectorXd Eval(const lf::mesh::Entity& entity);

 private:
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> _fe_space_p;
  Eigen::VectorXd _mu;
};

class progress_bar {
  static const auto overhead = sizeof " [100%]";
  std::ostream& os;
  const std::size_t bar_width;
  std::string message;
  const std::string full_bar;

 public:
  progress_bar(std::ostream& os, std::size_t line_width, std::string message_,
               const char symbol = '.')
      : os{os},
        bar_width{line_width - overhead},
        message{std::move(message_)},
        full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')} {
    if (message.size() + 1 >= bar_width ||
        message.find('\n') != std::string::npos) {
      os << message << '\n';
      message.clear();
    } else {
      message += ' ';
    }
    write(0.0);
  }

  progress_bar(const progress_bar&) = delete;
  progress_bar& operator=(const progress_bar&) = delete;

  ~progress_bar() {
    write(1.0);
    os << '\n';
  }

  void write(double fraction);
};

/* LIBRARY FUNCTIONS */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const lf::mesh::Entity& entity);

Eigen::VectorXd computeLumpedProjection(
    const lf::assemble::DofHandler& scal_dofh, const Eigen::VectorXd& mu,
    const lf::assemble::DofHandler& vec_dofh);

double computeL2Deviation(const lf::assemble::DofHandler& scal_dofh,
                          const Eigen::VectorXd& eta,
                          const lf::assemble::DofHandler& vec_dofh,
                          const Eigen::VectorXd& gamma);

Eigen::VectorXd solveBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p);

Eigen::VectorXd solveGradVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    const Eigen::VectorXd& mu, const lf::assemble::DofHandler& vec_dofh);

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p);

template <typename FUNCTOR_U>
Eigen::VectorXd interpolateData(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNCTOR_U&& u) {
  // Generate Lehrfem++ mesh functions out of the functors
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};

  Eigen::VectorXd dof_vector_u = lf::fe::NodalProjection(*fe_space_p, mf_u);

  return dof_vector_u;
};

}  // namespace ZienkiewiczZhuEstimator

#endif
