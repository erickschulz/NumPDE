/**
 * @file advectionsupg.cc
 * @brief NPDE homework AdvectionSUPG code
 * @author R. Hiptmair
 * @date July 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include "advectionsupg.h"

#include <lf/base/lf_assert.h>
#include <lf/base/ref_el.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/lagr_fe.h>
#include <math.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lf/mesh/test_utils/test_meshes.h"

namespace AdvectionSUPG {

/* SAM_LISTING_BEGIN_9 */
void enforce_boundary_conditions(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space,
    lf::assemble::COOMatrix<double>& A, Eigen::VectorXd& b) {
  const lf::mesh::utils::MeshFunctionGlobal mf_g_Gamma_in{
      [](const Eigen::Vector2d& x) { return std::pow(sin(M_PI * x(0)), 2); }};

  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(fe_space->Mesh(), false);
  // Loop over all edges
  for (const auto& edge : fe_space->Mesh()->Entities(1)) {
    LF_ASSERT_MSG(edge->RefEl() == lf::base::RefEl::kSegment(),
                  "Entity should be an edge");
    const lf::geometry::Geometry* geo_ptr = edge->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
    // Check if the edge lies on $\Gamma_{\mathrm{in}}$  (geometric test)
    if ((corners(0, 0) + corners(0, 1)) / 2. > 1. - 1e-5 ||
        (corners(1, 0) + corners(1, 1)) / 2. < 1e-5) {
      // Add the edge to the flagged entities
      bd_flags(*edge) = true;
    }
  }
  // Loop over all Points
  for (const auto& point : fe_space->Mesh()->Entities(2)) {
    LF_ASSERT_MSG(point->RefEl() == lf::base::RefEl::kPoint(),
                  "Entity should be an edge");
    const lf::geometry::Geometry* geo_ptr = point->Geometry();
    const Eigen::VectorXd coords = lf::geometry::Corners(*geo_ptr);
    // Check if the node lies on  $\Gamma_{\mathrm{in}}$ (geometric test)
    if (coords(0) > 1. - 1e-5 || coords(1) < 1e-5) {
      // Add the edge to the flagged entities
      bd_flags(*point) = true;
    }
  }
  auto flag_values{lf::fe::InitEssentialConditionFromFunction(
      *fe_space, bd_flags, mf_g_Gamma_in)};

  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&flag_values](lf::assemble::glb_idx_t dof_idx) {
        return flag_values[dof_idx];
      },
      A, b);
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveRotSUPG(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space) {
  // Velocity field
  auto vf = [](Eigen::Vector2d x) -> Eigen::Vector2d { return {-x[1], x[0]}; };
  lf::mesh::utils::MeshFunctionGlobal mf_v(vf);
  // Initialize the Matrix provider
  SUAdvectionElemMatrixProvider elmat_builder(mf_v);
  // Build the lhs Matrix
  lf::assemble::COOMatrix<double> A(fe_space->LocGlobMap().NumDofs(),
                                    fe_space->LocGlobMap().NumDofs());
  lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(),
                                      fe_space->LocGlobMap(), elmat_builder, A);
  Eigen::VectorXd b(fe_space->LocGlobMap().NumDofs());
  b.setZero();
  // Modify linear system to account for (essential) Dirichlet boundasry
  // conditions
  enforce_boundary_conditions(fe_space, A, b);
  // Solve the linear system of equations
  const Eigen::SparseMatrix<double> A_sparse = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  return solver.solve(b);
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void cvgL2SUPG() {
  // Generate triangular mesh of the unit square
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  builder.setBottomLeftCorner(Eigen::Vector2d{0., 0.})
      .setTopRightCorner(Eigen::Vector2d{1., 1.});

  // Initialize a vector that will hold all the approximated solutions
  std::vector<Eigen::VectorXd> u_h(8);
  // Initialize a vector that will hold all the approximate solutions
  std::vector<double> err(8);
  // Initialize a mesh function that will represent the analytical solution
  const lf::mesh::utils::MeshFunctionGlobal mf_u{
      [](const Eigen::Vector2d& x) -> double {
        if (x.norm() > 1.) return 0.;
        return std::pow(sin(M_PI * x.norm()), 2);
      }};

  std::cout << "N" << std::setw(8) << "|"
            << "L2 Error" << std::endl;

  // Compute error on 4 meshes of dyadically increasing resolution
  for (int i = 0; i < 4; i++) {
    // Refine mesh
    builder.setNumXCells(10 * pow(2, i)).setNumYCells(10 * pow(2, i));
    // Get the mesh
    auto mesh_p = builder.Build();
    // Set up finite element space
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
    // Compute the solution
    u_h[i] = AdvectionSUPG::solveRotSUPG(fe_space);

    const lf::fe::MeshFunctionFE mf_sol(fe_space, u_h[i]);
    // Compute the error
    err[i] = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 4));
    // Print the error in console
    std::cout << fe_space->LocGlobMap().NumDofs() << std::setw(20) << "|"
              << err[i] << std::endl;
    if (i == 3) {
      visSolution(fe_space, u_h[i]);
    }
  }
};
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_8 */
void visSolution(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space,
    Eigen::VectorXd& u) {
  const lf::fe::MeshFunctionFE mf_sol(fe_space, u);
  lf::io::VtkWriter vtk_writer(fe_space->Mesh(), "./solution.vtk");
  vtk_writer.WritePointData("solution", mf_sol);
}
/* SAM_LISTING_END_8 */

}  // namespace AdvectionSUPG
