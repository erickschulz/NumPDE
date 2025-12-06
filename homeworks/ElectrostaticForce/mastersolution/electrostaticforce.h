/**
 * @file electrostaticforce.h
 * @brief ElectrostaticForce
 * @author Erick Schulz
 * @date 27.11.2019
 * @copyright Developed at ETH Zurich
 */

#include <math.h>

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
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ElectrostaticForce {

/** @brief Compute the electrostatic force
 *       F(u) = (1/2) * int_{\Gamma_1} <grad(u),n> * grad(u) dS
 * from the potential directly using the trapezoidal rule. Recall that
 * for a partition x_0 < x_1 < ... < x_N of a curve C. The trapezoidal rule
 * reads which approximates
 *                         int_C f(x) dS
 * reads
 *      ( dx/2 ) * ( f(x_0) + 2*f(x_1) + ... + 2*f(x_N-1) + f(X_N) ).
 * Hence, when the curve is closed, that is when the endpoints x_0 = x_N are
 * the same, the quadrature rule simply reads
 *                      dx * sum_{i=0}^N f(x_i).
 * @return two dimensional force vector */
Eigen::Vector2d computeExactForce();

Eigen::VectorXd solvePoissonBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p);

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p);

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const lf::mesh::Entity& entity);

Eigen::Vector2d computeForceBoundaryFunctional(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::VectorXd approx_sol);

Eigen::Vector2d computeForceDomainFunctional(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::VectorXd approx_sol);

}  // namespace ElectrostaticForce
