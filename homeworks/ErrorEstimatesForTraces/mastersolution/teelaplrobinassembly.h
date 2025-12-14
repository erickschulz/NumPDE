/** @file
 * @brief NPDE ErrorEstimatesForTraces
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace ErrorEstimatesForTraces {

// Simplification of lenghty type names
using coord_t = Eigen::Vector2d;
using quad_rule_collection_t = std::map<lf::base::RefEl, lf::quad::QuadRule>;
using linear_lagrange = lf::uscalfe::FeSpaceLagrangeO1<double>;

// Function solving the Robin elliptic boundary value problem
Eigen::VectorXd solveBVP(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>&);

// Function for evaluating the integral of a function on the boundary
double bdFunctionalEval(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>&, Eigen::VectorXd&);

}  // namespace ErrorEstimatesForTraces
