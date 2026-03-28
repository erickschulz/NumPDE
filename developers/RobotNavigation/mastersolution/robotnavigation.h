/**
 * @file robotnavigation.h
 * @brief NPDE homework RobotNavigation code
 * @author Erick Schulz
 * @date 28.03.2026
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

namespace RobotNavigation {

using coord_t = Eigen::Vector2d;

/**
 * @brief Solve the Poisson BVP
 *
 *    -Delta(u) = 1 in room
 *            u = 0 at doors (Dirichlet)
 *        du/dn = 0 at walls (Neumann)
 *
 * on the mesh.
 */
Eigen::VectorXd solvePoissonBVP(const lf::io::GmshReader& reader);

/**
 * @brief Compute gradients of barycentric coordinate functions.
 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const lf::mesh::Entity& entity);

/**
 * @brief Element matrix provider for the vector FE mass matrix.
 *
 * Returns 6x6 (triangle) element mass matrices for the product space
 *
 *    V_h = S^0_1(M) x S^0_1(M),
 *
 * with DOF ordering
 *
 *    [gradb_x1, gradb_y1, gradb_x2, gradb_y2, gradb_x3, gradb_y3]
 *
 * per cell.
 */
class GradientProjectionMassMatrixProvider {
 public:
  explicit GradientProjectionMassMatrixProvider() = default;
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  Eigen::MatrixXd Eval(const lf::mesh::Entity& entity);
};

/**
 * @brief Element vector provider for the gradient projection RHS.
 *
 * Computes the element vector whose i-th entry is
 *
 *    ∫_K grad(u_h) · b_i dx
 *
 * where b_i are the local basis functions of V_h.
 * The integrals are approximated with the trapezoidal rule.
 *
 * @param u_h Coefficient vector of the scalar FE solution.
 */
class GradientProjectionRhsVectorProvider {
 public:
  explicit GradientProjectionRhsVectorProvider(
      const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
      const Eigen::VectorXd& u_h)
      : u_h_(u_h), fe_space_p_(fe_space_p) {}
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  Eigen::VectorXd Eval(const lf::mesh::Entity& entity);

 private:
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p_;
  Eigen::VectorXd u_h_;
};

/**
 * @brief Solve the gradient variational problem (3.9.2).
 *
 * Finds g_h in V_h such that
 *
 *    ∫ g_h · v_h = ∫ grad(u_h) · v_h
 *
 * for all v_h in V_h.
 *
 * Returns the V_h coefficient vector.
 */
Eigen::VectorXd solveGradientProjection(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    const Eigen::VectorXd& u_h, const lf::assemble::DofHandler& vec_dofh);

/**
 * @brief Trace the robot path from a starting point using RK4.
 *
 * Solves
 *
 *    dx/dt = -grad(u)(x(t))
 *
 * numerically with the 4th-order Runge-Kutta RK4.
 * Stops when u(x) < u_door_threshold (near a door) or after max_steps.
 *
 * @param mesh_p Pointer to the mesh.
 * @param dofh Scalar DOF handler.
 * @param sol Scalar FE solution coefficient vector.
 * @param vec_dofh Vector DOF handler (2 DOFs per node).
 * @param grad Gradient projection coefficient vector.
 * @param x_start Starting position of the robot.
 * @param dt Time step size.
 * @param max_steps Maximum number of steps.
 * @param u_door_threshold Stop when u(x) drops below this (near a door).
 * @return Nx2 matrix of path positions.
 */
Eigen::MatrixXd tracePath(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
                          const lf::assemble::DofHandler& dofh,
                          const Eigen::VectorXd& sol,
                          const lf::assemble::DofHandler& vec_dofh,
                          const Eigen::VectorXd& grad,
                          const Eigen::Vector2d& x_start, double dt,
                          int max_steps, double u_door_threshold);

}  // namespace RobotNavigation
