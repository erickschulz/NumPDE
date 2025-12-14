/**
 * @file semilinearellipticbvp.cc
 * @brief NPDE homework SemilinearEllipticBVP code
 * @author R. Hiptmair
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "semilinearellipticbvp.h"

#include <lf/fe/mesh_function_fe.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <iomanip>

namespace semilinearellipticbvp {
/* SAM_LISTING_BEGIN_3 */
void fixedPointNextIt(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    Eigen::VectorXd& mu_vec, Eigen::VectorXd& rhs_vec) {
  // Set up mesh functions for diffusion coefficient and reaction coefficient
  lf::mesh::utils::MeshFunctionGlobal mf_one(
      [](Eigen::Vector2d /*x*/) -> double { return 1.0; });
  const lf::fe::MeshFunctionFE mf_uh_prev(fes_p, mu_vec);
  const FunctionMFWrapper mf_coeff(mf_uh_prev, [](double xi) -> double {
    return (std::abs(xi) < 1.0E-16) ? 1.0 : std::sinh(xi) / xi;
  });

  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  const std::size_t N_dofs(dofh.NumDofs());
  // Assemble Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_one),
                                                      decltype(mf_coeff)>
      elmat_provider(fes_p, mf_one, mf_coeff);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);

  // Enforce homogeneous Dirichlet boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Elimination of degrees of freedom on the boundary. Also sets the
  // corresponding entries of rhs\_vec to zero.
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, rhs_vec);
  // Solve linear system
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  mu_vec = solver.solve(rhs_vec);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
}
/* SAM_LISTING_END_3 */

void testSolverSemilinearBVP(unsigned int reflevels) {
  // Manufactured solution of $-\Delta u + \sinh(u) = f$ with homogeneous
  // Dirichlet boundary conditions
  auto u = [](Eigen::Vector2d x) -> double {
    return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return M_PI *
           Eigen::Vector2d(std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]),
                           std::sin(M_PI * x[0]) * std::cos(M_PI * x[1]));
  };
  auto f = [&u](Eigen::Vector2d x) -> double {
    return 2.0 * M_PI * M_PI * u(x) + std::sinh(u(x));
  };
  // Wrap lambdas into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Triangular mesh hierarachy of unit square for testing
  // Adapted from Lehrfem++ \cppfile{homDir_linfe_demo.cc}
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << "\t Sequence of nested meshes used in test routine\n";
  multi_mesh.PrintInfo(std::cout);
  std::size_t L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeping error norms
  std::vector<std::tuple<int, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Compute finite-element solution of boundary value problem
    const Eigen::VectorXd sol_vec = solveSemilinearBVP(fe_space, f);

    // Compute error norms
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(mesh_p->NumEntities(2), L2err, H1serr);
  }
  // Output table of errors
  std::cout << "\t Table of error norms" << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 error" << std::setw(16) << "H1 error" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (const auto& err : errs) {
    auto [N, l2err, h1serr] = err;
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::setw(16) << h1serr << std::endl;
  }
}

/* SAM_LISTING_BEGIN_4 */
void newtonNextIt(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    Eigen::VectorXd& mu_vec, Eigen::VectorXd& rhs_vec) {
  // We need x->0, x->1, and x->uh(x) as MeshFunctions
  lf::mesh::utils::MeshFunctionConstant mf_one(1.);
  lf::mesh::utils::MeshFunctionConstant mf_zero(0.);
  const lf::fe::MeshFunctionFE mf_uh_prev(fes_p, mu_vec);

  // We need x->cosh(uh\_prev(x)) and x->-sinh(uh\_prev(x)) as MeshFunctions
  const FunctionMFWrapper mf_cosh_uh_prev(
      mf_uh_prev, [](double xi) -> double { return std::cosh(xi); });
  const FunctionMFWrapper mf_sinh_uh_prev(
      mf_uh_prev, [](double xi) -> double { return -std::sinh(xi); });

  // We define some variables for easy access
  const lf::mesh::Mesh& mesh{*(fes_p->Mesh())};
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  const std::size_t N_dofs(dofh.NumDofs());

  // Assemble matrix A, representing (grad uh, grad vh) + (cosh(uh\_prev)*uh,
  // vh) in COO format
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_one),
                                                      decltype(mf_cosh_uh_prev)>
      elmat_provider(fes_p, mf_one, mf_cosh_uh_prev);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);

  // Assemble matrix M\_crs, representing (grad uh, grad vh) in
  // Eigen::SparseMatrix format
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_one),
                                                      decltype(mf_zero)>
      elmat_provider_M(fes_p, mf_one, mf_zero);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider_M, M);
  Eigen::SparseMatrix<double> M_crs = M.makeSparse();

  // Construct the right-hand side of the system
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();

  // Add (f,vh)
  phi += rhs_vec;

  // Add (-sinh(uh\_prev),vh)
  lf::uscalfe::ScalarLoadElementVectorProvider<double,
                                               decltype(mf_sinh_uh_prev)>
      elvec_builder(fes_p, mf_sinh_uh_prev);
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Add -(grad uh\_prev, grad vh)
  phi -= M_crs * mu_vec;

  // Enforce homogeneous Dirichlet boundary conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Elimination of degrees of freedom on the boundary. Also sets the
  // corresponding entries of rhs\_vec to zero.
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);

  // Convert matrix A in COO format to A\_crs in Eigen::SparseMatrix format
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve the system and add the solution to mu\_vec
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  mu_vec += solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
};
/* SAM_LISTING_END_4 */

}  // namespace semilinearellipticbvp
