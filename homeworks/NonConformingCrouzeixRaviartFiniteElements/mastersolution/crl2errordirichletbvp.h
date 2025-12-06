/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   18.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
#define NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H

#include <lf/io/io.h>
#include <lf/mesh/mesh.h>

#include <cmath>
#include <string>

#include "crdirichletbvp.h"
#include "crfespace.h"
#include "crl2error.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
double L2errorCRDiscretizationDirichletBVP(const std::string& filename) {
  double l2_error;

// TODO: task 2-14.x)
  // Right-hand-side source function
  auto f = [](Eigen::Vector2d x) -> double {
    return (2. * M_PI * M_PI + x.prod()) * std::sin(M_PI * x(0)) *
           std::sin(M_PI * x(1));
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return x.prod(); };
  // Analytic solution
  auto u = [](Eigen::Vector2d x) -> double {
    return std::sin(M_PI * x(0)) * std::sin(M_PI * x(1));
  };

  // Read mesh from file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), filename);

  // Build CR FE space
  auto fe_space = std::make_shared<CRFeSpace>(reader.mesh());

  // Solve homogeneous Dirichlet problem
  Eigen::VectorXd mu = solveCRDirichletBVP(fe_space, gamma, f);
  // Compute L2 norm of error
  l2_error = computeCRL2Error(fe_space, mu, u);
  return l2_error;
}
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
