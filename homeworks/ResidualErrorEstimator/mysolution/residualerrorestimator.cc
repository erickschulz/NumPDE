/**
 * @file
 * @brief NPDE homework ResidualErrorEstimator
 * @author Ralf Hiptmair
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include "residualerrorestimator.h"

#include <utility>

namespace REE {

/* SAM_LISTING_BEGIN_2 */
dataDiscreteBVP::dataDiscreteBVP(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                                 std::function<double(Eigen::Vector2d)> alpha,
                                 std::function<double(Eigen::Vector2d)> f)
    : pwlinfespace_p_(
          std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p)),
      mf_f_(f),
      mf_alpha_(mesh_p, alpha) {
  std::cout << "Discrete BVP built" << std::endl;
}
/* SAM_LISTING_END_2 */

Eigen::VectorXd solveBVP(const dataDiscreteBVP& disc_bvp) {
  // For conveneicne we set up references to essential objects for FE
  // discretization in the lowest-order Lagrangian finite element space
  const lf::uscalfe::FeSpaceLagrangeO1<double>& linfespc{
      *disc_bvp.pwlinfespace_p_};
  // The underlying finite-element mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{linfespc.Mesh()};
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{linfespc.LocGlobMap()};
  // Dimension of finite element space, number of unknowns
  const lf::base::size_type N_dofs(dofh.NumDofs());

  // I: Assembly of full Galerkin matrix
  // Object for sparse matrix to be filled by cell-oriented assembly
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Provider object for element matrices for scalar linear second-order pure
  // diffusion operator with variable diffusion coefficient. Uses numerical
  // quadrature of order 3 and, thus, computes exact element matrices for
  // locally constant diffusion coefficient.
  lf::fe::DiffusionElementMatrixProvider<double, decltype(disc_bvp.mf_alpha_)>
      elmat_builder(disc_bvp.pwlinfespace_p_, disc_bvp.mf_alpha_);
  // Invoke cell-oriented assembly of the finite-element Galerkin matrix
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // II: Assembly of right-hand-side vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Assemble volume part of right-hand side vector depending on the source
  // function f.
  // Initialize object taking care of local computations on all cells.
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(disc_bvp.mf_f_)>
      elvec_builder(disc_bvp.pwlinfespace_p_, disc_bvp.mf_f_);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // III: Enforce zero Dirichlet (essential) boundary conditions
  // Create a predicate selecting nodes on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t dof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(dof_idx)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  // Examine return status of solver in case the matrix is singular
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  return sol_vec;
}

/* SAM_LISTING_BEGIN_3 */
lf::mesh::utils::CodimMeshDataSet<double> volumeResiduals(
    const dataDiscreteBVP& disc_bvp, const Eigen::VectorXd& /*u_vec*/) {
  // Get pointer to underlying mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      disc_bvp.pwlinfespace_p_->Mesh();
  // Mesh data set for returning the volume residuals
  lf::mesh::utils::CodimMeshDataSet<double> vol_res(mesh_p, 0, 0.0);

  //====================
  // Your code goes here
  //====================
  return vol_res;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
lf::mesh::utils::CodimMeshDataSet<double> edgeResiduals(
    const dataDiscreteBVP& disc_bvp, const Eigen::VectorXd& u_vec) {
  // Get pointer to underlying mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      disc_bvp.pwlinfespace_p_->Mesh();
  // Mesh data set for returning the volume residuals
  lf::mesh::utils::CodimMeshDataSet<double> edge_res(mesh_p, 1, 0.0);
  //====================
  // Your code goes here
  //====================
  return edge_res;
}
/* SAM_LISTING_END_4 */

std::tuple<double, double, double> solveAndEstimate(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Note: the mesh must cover the unit square for this test setting !

  // Define homogeneous Dirichlet boundary value problem
  // Diffusion coefficient function
  std::function<double(Eigen::Vector2d)> alpha =
      [](Eigen::Vector2d x) -> double { return (1.0 /* + x.squaredNorm() */); };
  /*
  // Right-hand side source function
  std::function<double(Eigen::Vector2d)> f = [](Eigen::Vector2d x) -> double {
    return (std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]));
  };
  // For testing: exact solution
  std::function<double(Eigen::Vector2d)> u_exact =
      [&f](Eigen::Vector2d x) -> double { return (f(x) / (M_PI * M_PI)); };
  // For testing: gradient of exact solution
  std::function<Eigen::Vector2d(Eigen::Vector2d)> grad_u =
      [](Eigen::Vector2d x) -> Eigen::Vector2d {
    const double den = M_PI;
    return ((Eigen::Vector2d() << std::cos(M_PI * x(0)) * std::sin(M_PI * x(1)),
             std::sin(M_PI * x(0)) * std::cos(M_PI * x(1)))
                .finished()) /
           den;
  };
  */
  auto u_exact = [](Eigen::Vector2d x) -> double {
    return (x[0] * (1.0 - x[0]) * x[1] * (1.0 - x[1]));
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return ((Eigen::Vector2d() << (1 - 2 * x[0]) * x[1] * (1 - x[1]),
             x[0] * (1 - x[0]) * (1 - 2 * x[1]))
                .finished());
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return 2 * (x[1] * (1 - x[1]) + x[0] * (1 - x[0]));
  };

  // Lambdas have to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u_exact};
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};

  // Defining the discretized boundary value problem including the
  // finite-element space
  const dataDiscreteBVP disc_bvp(mesh_p, alpha, f);

  // Compute basis expansion coefficient vector of finite-element solution
  Eigen::VectorXd mu{solveBVP(disc_bvp)};

  // Compute error norms
  // create *mesh_p functions representing solution / gradient of solution
  const lf::fe::MeshFunctionFE mf_sol(disc_bvp.pwlinfespace_p_, mu);
  const lf::fe::MeshFunctionGradFE mf_grad_sol(disc_bvp.pwlinfespace_p_, mu);
  // compute errors with 3rd order quadrature rules, which is sufficient for
  // piecewise linear finite elements
  double L2err =  // NOLINT
      std::sqrt(lf::fe::IntegrateMeshFunction(
          *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
  double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
      *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));

  std::cout << "Mesh (" << mesh_p->NumEntities(0) << " cells, "
            << mesh_p->NumEntities(1) << " edges, " << mesh_p->NumEntities(2)
            << " nodes): L2err = " << L2err << ", H1serr = " << H1serr
            << std::endl;

  // Evaluate a-posteriori error estimator
  // Compute volume contributions
  lf::mesh::utils::CodimMeshDataSet<double> vol_res{
      volumeResiduals(disc_bvp, mu)};
  // Compute edge contributions
  lf::mesh::utils::CodimMeshDataSet<double> ed_res{edgeResiduals(disc_bvp, mu)};

  // Sum volume residuals and edge residuals
  double eta_vol = 0.0;
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    eta_vol += vol_res(*cell);
  }
  double eta_ed = 0.0;
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    eta_ed += ed_res(*edge);
  }

  std::cout << "Estimated error = " << eta_vol + eta_ed << ", vol = " << eta_vol
            << ", edge = " << eta_ed << std::endl;

  return {L2err, H1serr, eta_vol + eta_ed};
}  // end solveAndEstimate

}  // namespace REE
