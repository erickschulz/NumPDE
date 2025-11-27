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

#if SOLUTION
  // Edge midpoint quadrature rule for triangles, which is exact for quadratic
  // polynomials and allows the exact computation of the volume residual
  // contributions for piecewise linear source functions.
  const lf::quad::QuadRule qr{lf::quad::make_TriaQR_EdgeMidpointRule()};
  const Eigen::VectorXd qw{qr.Weights()};  // quadrature weigth vector

  // Run over all cells of the mesh
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Implemented for triangles only");
    // Obtain information about the shape of the cell
    const lf::geometry::Geometry& geo{*(cell->Geometry())};
    // Ddetermine size of triangle (length of longest edge)
    const Eigen::MatrixXd corners{lf::geometry::Corners(geo)};
    const double h0 = (corners.col(1) - corners.col(0)).norm();
    const double h1 = (corners.col(2) - corners.col(1)).norm();
    const double h2 = (corners.col(0) - corners.col(2)).norm();
    const double h_K = std::max({h0, h1, h2});

    // I: Locally integrate the square of the source function
    // Gramian determinants
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(qr.Points())};
    // Fetch values of source function in quadrature nodes. Note that local
    // reference coordinates have to be passed to the evaluation operator of the
    // MeshFunction.
    const std::vector<double> f_vals{disc_bvp.mf_f_(*cell, qr.Points())};
    // Evaluation of quadrature formula
    double f_sq_int = 0.0;
    for (lf::base::size_type l = 0; l < qr.NumPoints(); ++l) {
      f_sq_int += qw[l] * f_vals[l] * f_vals[l] * gram_dets[l];
    }
    // II: Fetch prefactor depending on diffusion coefficient
    // We need a "dummy point" in the reference element in order to call the
    // evaluation operator for the MeshFunction.
    const Eigen::MatrixXd dummy = Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0);
    const double alpha_K = disc_bvp.mf_alpha_(*cell, dummy)[0];
    LF_ASSERT_MSG(alpha_K > 0,
                  "Diffusion coefficients must be strictly positive!");
    const double mu_K = h_K * h_K / alpha_K;
    // Store cell contribution
    vol_res(*cell) = f_sq_int * mu_K;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
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
#if SOLUTION
  // Step I: Create an auxiliary MeshDataSet storing non-normalized (!) edge
  // normals = edge direction vectors turned by 90 degrees
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> edge_normals(mesh_p, 1);
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> edge_startpt(mesh_p, 1);
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    // Obtain information about the shape of the edge
    const lf::geometry::Geometry& geo{*(edge->Geometry())};
    const Eigen::MatrixXd corners{lf::geometry::Corners(geo)};
    // Starting point of the edge
    edge_startpt(*edge) = corners.col(0);
    // Direction vector of the edge
    const Eigen::Vector2d dir = corners.col(1) - corners.col(0);
    // Rotate counterclockwise by 90 degrees an d store
    edge_normals(*edge) = Eigen::Vector2d(dir(1), -dir(0));
  }
  // Step II: Compute the gradient of the finite element solution, which
  // LehrFEM++ can do for us.
  const lf::fe::MeshFunctionGradFE grad_uh(disc_bvp.pwlinfespace_p_, u_vec);

  // Step III: Visit the triangles of the mesh, fetch the constant(!) gradient
  // of the finite-element solution, multiply it with the edge normals and the
  // coefficient and sum the results for every edge not loacted on the
  // bounmdary. This summation is done in an auxiliary MeshDataSet.
  lf::mesh::utils::CodimMeshDataSet<double> alpha_max(mesh_p, 1, 0.0);
  lf::mesh::utils::CodimMeshDataSet<double> edge_flux_jump(mesh_p, 1, 0.0);
  lf::mesh::utils::CodimMeshDataSet<bool> bd_ed_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Implemented for triangles only");
    // Retrieve constant diffusion coefficient
    const Eigen::MatrixXd dummy = Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0);
    const double alpha_K = disc_bvp.mf_alpha_(*cell, dummy)[0];
    LF_ASSERT_MSG(alpha_K > 0,
                  "Diffusion coefficients must be strictly positive!");
    // Fetch scaled gradient of the finite element solution
    const Eigen::Vector2d nablau_K = (grad_uh(*cell, dummy)[0]) * alpha_K;
    /* Debugging output
    std::cout << "cell " << mesh_p->Index(*cell)
              << ": alpha*grad u_h = " << nablau_K.transpose() << std::endl;
    */
    // Obtain shape of cell
    const lf::geometry::Geometry& cell_geo{*(cell->Geometry())};
    const Eigen::MatrixXd cell_vert{lf::geometry::Corners(cell_geo)};
    // Visit all three edges of the triangle
    std::span<const lf::mesh::Entity* const> edges{cell->SubEntities(1)};
    LF_ASSERT_MSG(edges.size() == 3, "Triangle must have three edges!");
    for (int l = 0; l < 3; ++l) {
      if (!bd_ed_flags(*edges[l])) {
        // Determine orientation of edge normal
        const lf::mesh::Entity& edge = *edges[l];
        const Eigen::Vector2d& normal{edge_normals(edge)};
        const int ori =
            (normal.dot(cell_vert.col((l + 2) % 3) - edge_startpt(edge)) > 0)
                ? 1
                : -1;
        edge_flux_jump(edge) += (ori * nablau_K.dot(normal));
        /* Debugging output
        std::cout << "edge " << mesh_p->Index(edge) << ", r.ori = " << ori
                  << ": edge jump = " << edge_flux_jump(edge) << std::endl;
        */
      }
      if (alpha_K > alpha_max(*edges[l])) {
        alpha_max(*edges[l]) = alpha_K;
      }
    }
  }
  // Now we have all required information in the auxiliary MeshDataSets
  // and we traverse the edges again and compute the scaled jump norms.
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    const double ed_flux = edge_flux_jump(*edge);
    edge_res(*edge) = (ed_flux * ed_flux) / alpha_max(*edge);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
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
