/**
 * @file electrostaticforce.cc
 * @brief ElectrostaticForce
 * @author Erick Schulz
 * @date 27.11.2019
 * @copyright Developed at ETH Zurich
 */

#include "electrostaticforce.h"

namespace ElectrostaticForce {

Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const lf::mesh::Entity& entity) {
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(entity.Geometry()));

  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();

  return X.inverse().block<2, 3>(1, 0);
}  // gradbarycoordinates

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d computeExactForce() {
  Eigen::Vector2d force;  // return vector
  // Real overkill quadrature!
  unsigned int N = 1e6;  // nb. of quadrature points
  // Data
  double r = 4.0 / 15.0;  // radius of the circle
  Eigen::Vector2d a(-16.0 / 15.0, 0.0);
  Eigen::Vector2d b(-1.0 / 15.0, 0.0);
  // Auxiliary variables
  Eigen::Vector2d x;     // euclidean coordinates
  Eigen::Vector2d n;     // normal vector at the coordinates
  Eigen::Vector2d grad;  // gradient of the exact solution
  double angle;          // interior angle

  // Compute the force using the trapezoidal rule
  force.setZero();  // initially zero
  for (unsigned int i = 0; i < N; i++) {
    // I. Find the euclidean coordinates of integration node
    angle = i * (2 * M_PI) / N;
    x(0) = r * std::cos(angle);
    x(1) = r * std::sin(angle);
    // II. Evaluate the gradient at the coordinates
    grad = ((x - a) / (x - a).squaredNorm() - (x - b) / (x - b).squaredNorm()) /
           std::log(2.0);
    // III. Evaluate the normal vector at the coordinates
    n = -x / x.norm();
    // IV. Evaluate full integrand
    force += grad.dot(n) * grad;
  }
  // V. Scale the result of the summation
  force *= r * M_PI / N;
  return force;
}  // end computeExactForce
/* SAM_LISTING_END_1 */

Eigen::VectorXd solvePoissonBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p) {
  Eigen::VectorXd approx_sol;  // to return

  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Producing Dirichlet data
  /* SAM_LISTING_BEGIN_2 */
  auto mf_bd_values =
      lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
        if (x.norm() < 0.27) {
          return 1.0;
        }
        return 0.0;
      });
  /* SAM_LISTING_END_2 */
  /* SAM_LISTING_BEGIN_3 */
  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Computing volume matrix for negative Laplace operator
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // II. IMPOSING ESSENTIAL BOUNDARY CONDITIONS

  // Obtain arrays of boolean flags for the edges and nodes of the mesh, 'true'
  // indicates that the edge or node lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Determine the fixed dofs on the boundary and their values
  auto flag_values{lf::fe::InitEssentialConditionFromFunction(
      *fe_space_p, bd_flags, mf_bd_values)};
  // II.ii Eliminate Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&flag_values](lf::assemble::glb_idx_t dof_idx) {
        return flag_values[dof_idx];
      },
      A, phi);
  /* SAM_LISTING_END_3 */

  // II : SOLVING  THE LINEAR SYSTEM
  //  Convert COO format to CRS using Eigen's internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();
  // II.i : Setting up Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // II.ii : Solving
  approx_sol = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  LF_VERIFY_MSG(approx_sol.size() == N_dofs, "Solution has wrong size");
  return approx_sol;
}  // solvePoissonBVP

/* SAM_LISTING_BEGIN_4 */
Eigen::Vector2d computeForceBoundaryFunctional(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::VectorXd approx_sol) {
  Eigen::Vector2d approx_force;  // to return

  // MESH AND FINITE ELEMENTS SPACE DATA
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Obtain arrays of boolean flags for the edges and nodes of the mesh, 'true'
  // indicates that the edge or node lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // This predicate returns 'true' if the edge belongs to the interior boundary
  auto interior_bd_flags = [&bd_flags](const lf::mesh::Entity& edge) -> bool {
    if (bd_flags(edge)) {
      auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
      if (endpoints.col(0).norm() < 0.27) {
        return true;
      }
    }
    return false;
  };

  // INTEGRATION TOOLS
  double edge_length;
  Eigen::VectorXd normal_vec;
  Eigen::Matrix2d rotation_mat;
  rotation_mat << 0, 1, -1, 0;             // rotates a 2d vec by 90 deg.
  Eigen::Matrix<double, 2, 3> elgrad_mat;  // gradients of basis elements
  Eigen::Vector2d loc_grad_approx_sol;  // gradient expansion coeff FE solution

  // PERFORMING INTEGRATION
  approx_force.setZero();
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    for (const lf::mesh::Entity* edge : cell->SubEntities(1)) {
      if (interior_bd_flags(*edge)) {
        auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
        edge_length = (endpoints.col(1) - endpoints.col(0)).norm();

        // I COMPUTING NORMAL VECTOR
        // I.i Finding a unit vector perpendicular to the edge
        normal_vec =
            rotation_mat * (endpoints.col(1) - endpoints.col(0)).normalized();
        // I.ii The normal vector must be pointing outward
        if (normal_vec.dot(endpoints.col(0)) > 0) {
          normal_vec *= -1.0;
        }

        // II COMPUTING GRADIENT OF THE FE APPROXIMATE SOLUTION
        // II.i Obtain global FE indices of the vertices
        const auto dof_idx_vec = dofh.GlobalDofIndices(*cell);
        // II.ii Obtain the gradients of the barycentric coordinate functions
        elgrad_mat = gradbarycoordinates(*cell);
        // II.iii Compute the gradient of the passed coefficient vector
        loc_grad_approx_sol = elgrad_mat.col(0) * approx_sol(dof_idx_vec[0]) +
                              elgrad_mat.col(1) * approx_sol(dof_idx_vec[1]) +
                              elgrad_mat.col(2) * approx_sol(dof_idx_vec[2]);

        // III SUMMING LOCAL CONTRIBUTION
        approx_force += edge_length * loc_grad_approx_sol.dot(normal_vec) *
                        loc_grad_approx_sol;
      }
    }
  }
  approx_force *= 0.5;
  return approx_force;
}  // computeForceFunctional
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::Vector2d computeForceDomainFunctional(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space_p,
    Eigen::VectorXd approx_sol) {
  Eigen::Vector2d approx_force;  // to return

  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};

  Eigen::Vector2d a(-16.0 / 15.0, 0.0);
  Eigen::Vector2d b(-1.0 / 15.0, 0.0);
  auto grad_uExact = [&a, &b](Eigen::VectorXd x) -> Eigen::Vector2d {
    return ((x - a) / (x - a).squaredNorm() - (x - b) / (x - b).squaredNorm()) /
           std::log(2.0);
  };

  // INTEGRATION TOOLS
  Eigen::MatrixXd stress_tensor(2, 2);     // maxwell stress tensor
  Eigen::MatrixXd mid_pts(2, 3);           // midpoints of the edges
  Eigen::Matrix<double, 2, 3> elgrad_mat;  // gradients of basis elements
  Eigen::Vector2d loc_grad_approx_sol;  // gradient expansion coeff FE solution
  Eigen::MatrixXd Id = Eigen::Matrix<double, 2, 2>::Identity();

  // PERFORMING INTEGRATION
  approx_force.setZero();
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    auto endpoints = lf::geometry::Corners(*(cell->Geometry()));
    mid_pts.col(0) = 0.5 * (endpoints.col(0) + endpoints.col(1));
    mid_pts.col(1) = 0.5 * (endpoints.col(1) + endpoints.col(2));
    mid_pts.col(2) = 0.5 * (endpoints.col(2) + endpoints.col(1));

    double area = lf::geometry::Volume(*(cell->Geometry()));

    // II COMPUTING GRADIENT OF THE FE APPROXIMATE SOLUTION
    // II.i Obtain global FE indices of the vertices
    const auto dof_idx_vec = dofh.GlobalDofIndices(*cell);
    // II.ii Obtain the gradients of the barycentric coordinate functions
    elgrad_mat = gradbarycoordinates(*cell);
    // II.iii Compute the gradient of the passed coefficient vector
    loc_grad_approx_sol = elgrad_mat.col(0) * approx_sol(dof_idx_vec[0]) +
                          elgrad_mat.col(1) * approx_sol(dof_idx_vec[1]) +
                          elgrad_mat.col(2) * approx_sol(dof_idx_vec[2]);

    // III ASSEMBLING MAXWELL STRESS TENSOR
    stress_tensor = loc_grad_approx_sol * loc_grad_approx_sol.transpose() -
                    0.5 * loc_grad_approx_sol.squaredNorm() * Id;

    // IV SUMMING LOCAL CONTRIBUTION
    approx_force += (area / 3.0) * stress_tensor *
                    (grad_uExact(mid_pts.col(0)) + grad_uExact(mid_pts.col(1)) +
                     grad_uExact(mid_pts.col(2)));
  }
  return approx_force;
}  // computeForceDomainFunctional
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p) {
  double mesh_size = 0.0;

  // Find maximal edge length
  double edge_length;
  // Loop over all edges of the mesh
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    // Compute the length of the edge
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }
  return mesh_size;
};  // getMeshSize
/* SAM_LISTING_END_6 */

}  // namespace ElectrostaticForce
