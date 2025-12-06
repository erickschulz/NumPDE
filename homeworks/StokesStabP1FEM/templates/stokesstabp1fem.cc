/**
 * @file stokesstabp1fem.cc
 * @brief NPDE homework StokesStabP1FEM code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#include "stokesstabp1fem.h"

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_global.h>

namespace StokesStabP1FEM {

/* SAM_LISTING_BEGIN_1 */
P1StabFEMElementMatrixProvider::ElemMat P1StabFEMElementMatrixProvider::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  LF_VERIFY_MSG(cell.Geometry()->isAffine(),
                "Triangle must have straight edges");
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Compute gradients of barycentric coordinate functions and store them in the
  // columns of the $2\times 3$-matrix G
  // clang-format off
  const Eigen::MatrixXd dpt = (Eigen::MatrixXd(2, 1) << 0.0, 0.0).finished();
  const Eigen::Matrix<double, 2, 3> G =
      cell.Geometry()->JacobianInverseGramian(dpt).block(0, 0, 2, 2) *
    (Eigen::Matrix<double, 2, 3>(2,3) << -1, 1, 0,
                                         -1, 0, 1).finished();
  // clang-format on
  // Compute the element matrix  for $-\Delta$ and $\cob{\Cs^0_1}$.
  // See also \lref{mc:ElementMatrixLaplLFE}.
  Eigen::Matrix<double, 3, 3> L = area * G.transpose() * G;

  // Do not forget to set all non-initialized entries to zero
  MK_.setZero();
  // Arrays for local index remapping
  const std::array<Eigen::Index, 3> vx_idx{0, 3, 6};
  const std::array<Eigen::Index, 3> vy_idx{1, 4, 7};
  const std::array<Eigen::Index, 3> p_idx{2, 5, 8};
  // Distribute the entries of L to the final element matrix
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      MK_(vx_idx[i], vx_idx[j]) = MK_(vy_idx[i], vy_idx[j]) = L(i, j);
    }
  }
  // Fill entries related to the B-blocks of the element matrix
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // \prbeqref{eq:BK} with 3-point edge midpoint quadrature rule!
      MK_(p_idx[i], vx_idx[j]) = MK_(vx_idx[j], p_idx[i]) =
          G(0, j) * area / 3.0;
      MK_(p_idx[i], vy_idx[j]) = MK_(vy_idx[j], p_idx[i]) =
          G(1, j) * area / 3.0;
    }
  }
  // Fill $3\times 3$ block due to stabilization terms
  // Element mass matrix for $\LFE$
  /* A more compact option for implementation
  const Eigen::Matrix3d SK =
   (area / 12.0) *
      (Eigen::Matrix3d() << 2, -1, -1, -1, 2, -1, -1, -1, 2).finished();
  */
  Eigen::Matrix3d SK;
  SK(0, 0) = (area / 12.0) * 2;
  SK(1, 1) = (area / 12.0) * 2;
  SK(2, 2) = (area / 12.0) * 2;
  SK(0, 1) = SK(1, 0) = -(area / 12.0);
  SK(0, 2) = SK(2, 0) = -(area / 12.0);
  SK(2, 1) = SK(1, 2) = -(area / 12.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      MK_(p_idx[i], p_idx[j]) = -SK(i, j);
    }
  }
  return MK_;
}
/* SAM_LISTING_END_1 */

lf::assemble::COOMatrix<double> buildP1StabFEMGalerkinMatrix(
    const lf::assemble::DofHandler& dofh) {
  // Total number of FE d.o.f.s without Lagrangian multiplier
  lf::assemble::size_type n = dofh.NumDofs();
  // Full Galerkin matrix in triplet format taking into account the zero mean
  // constraint on the pressure.
  lf::assemble::COOMatrix<double> A(n + 1, n + 1);
  // Set up computation of element matrix
  P1StabFEMElementMatrixProvider emp{};
  // Assemble \cor{full} Galerkin matrix for Taylor-Hood FEM
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, emp, A);

  // Add bottom row and right column corresponding to Lagrange multiplier
  // You cannot use AssembleMatrixLocally() because the DofHandler does
  // not know about this extra unknown.
  // Do cell-oriented assembly "manually"
  for (const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Only implemented for triangles");
    // Obtain area of triangle
    const double area = lf::geometry::Volume(*cell->Geometry());
    // The pressure GSFs are associated with the nodes
    const std::span<const lf::mesh::Entity* const> nodes{cell->SubEntities(2)};
    // Loop over nodes
    for (const lf::mesh::Entity* node : nodes) {
      // Area of the cell
      // Obtain index of tent function associated with node
      // All indices of global shape functions sitting at node
      std::span<const lf::assemble::gdof_idx_t> dof_idx{
          dofh.InteriorGlobalDofIndices(*node)};
      LF_ASSERT_MSG(dof_idx.size() == 3, "Node must carry 3 dofs!");
      // The index of the pressure global shape function is the third one
      const lf::assemble::gdof_idx_t tent_idx = dof_idx[2];
      A.AddToEntry(n, tent_idx, area / 3.0);
      A.AddToEntry(tent_idx, n, area / 3.0);
    }
  }
  // Rely on return value optimization
  return A;
}

void testCvgP1StabFEM(unsigned int refsteps) {
  using namespace std::numbers;
  // ********** Part I: Manufactured solution  **********
  // Analytic solution for velocity and pressure
  // Divergence-free velocity field with vanishing Dirichlet trace
  auto v_ex = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(
        -std::pow(std::sin(pi * x[0]), 2) * std::sin(2 * pi * x[1]),
        std::sin(2 * pi * x[0]) * std::pow(std::sin(pi * x[1]), 2));
  };
  // Pressure with vanishing mean
  auto p_ex = [](Eigen::Vector2d x) -> double { return (x[0] + x[1] - 1); };
  // Right-hand side forcing field
  auto f = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-2 * pi * pi * (1.0 - 2 * std::cos(2 * pi * x[0])) *
                                   std::sin(2 * pi * x[1]) -
                               1.0,
                           -2 * pi * pi * std::sin(2 * pi * x[0]) *
                                   (2 * std::cos(2 * pi * x[1]) - 1.0) -
                               1.0);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);

  auto grad_v1 = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(
        -pi * std::sin(2 * pi * x[0]) * std::sin(2 * pi * x[1]),
        -2 * pi * std::cos(2 * pi * x[1]) * std::pow(std::sin(pi * x[0]), 2));
  };
  auto grad_v2 = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(
        2 * pi * std::cos(2 * pi * x[0]) * std::pow(std::sin(pi * x[1]), 2),
        pi * std::sin(2 * pi * x[0]) * std::sin(2 * pi * x[1]));
  };
  // ********** Part II: Loop over sequence of meshes **********
  // Generate a small unstructured triangular mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);

  // Number of levels
  const int L = multi_mesh.NumLevels();

  // Table of various error norms
  std::vector<std::tuple<size_t, double, double, double, double, double>> errs;
  for (int level = 0; level < L; ++level) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Define Lagranggian FE spaces for piecewise linear approximation
    auto fes_o1_ptr =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(lev_mesh_p);
    // Fetch dof handler for the components of the velocity
    const lf::assemble::DofHandler& dofh_u = fes_o1_ptr->LocGlobMap();
    //  Fetch dof handler for the pressure (the same as for velocity)
    const lf::assemble::DofHandler& dofh_p = fes_o1_ptr->LocGlobMap();
    // ********** Part III: Solving on a single mesh **********
    // Initialize dof handler for P1 FEM for both velocity and pressure
    lf::assemble::UniformFEDofHandler dofh(lev_mesh_p,
                                           {{lf::base::RefEl::kPoint(), 3},
                                            {lf::base::RefEl::kSegment(), 0},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 0}});
    LF_ASSERT_MSG(dofh.NumDofs() == 2 * dofh_u.NumDofs() + dofh_p.NumDofs(),
                  "No dof mismatch");
    // Total number of d.o.f. in monolithic FE spaces
    size_t n = dofh.NumDofs();
    // Build and solve the linear system with trace of the exact velocity
    // solution as Dirichlet data.
    const Eigen::VectorXd res = solveP1StabFEMStokesBVP(dofh, mf_f);
    std::cout << "Computing with " << n << " d.o.f.s, solving .. "
              << std::flush;

    // Coefficient vectors for the first and second component of the velocity
    Eigen::VectorXd coeff_vec_u1 = Eigen::VectorXd::Zero(dofh_u.NumDofs());
    Eigen::VectorXd coeff_vec_u2 = Eigen::VectorXd::Zero(dofh_u.NumDofs());
    // Coefficient vector for the pressure
    Eigen::VectorXd coeff_vec_p = Eigen::VectorXd::Zero(dofh_p.NumDofs());

    // ********** Part IV: Compute error norms **********
    // Remapping dofs in order to be able to use MeshFunctionFE
    for (auto e : lev_mesh_p->Entities(2)) {
      // Global indices for u1, u2 for the respective vertex or edge
      auto glob_idxs = dofh.InteriorGlobalDofIndices(*e);
      auto glob_idx_o1 = dofh_u.InteriorGlobalDofIndices(*e)[0];
      // Extract the correct elements for the coefficient vector of the
      // components of u and the pressure p
      coeff_vec_u1[glob_idx_o1] = res[glob_idxs[0]];
      coeff_vec_u2[glob_idx_o1] = res[glob_idxs[1]];
      coeff_vec_p(glob_idx_o1) = res[glob_idxs[2]];
    }
    // Variables for storing the error norms
    double L2err_u1, L2err_u2, H1err_u1, H1err_u2, L2err_p;

    // Define finite-element mesh functions
    const lf::fe::MeshFunctionFE mf_o2_u1(fes_o1_ptr, coeff_vec_u1);
    const lf::fe::MeshFunctionFE mf_o2_u2(fes_o1_ptr, coeff_vec_u2);
    const lf::fe::MeshFunctionFE mf_o1_p(fes_o1_ptr, coeff_vec_p);
    const lf::fe::MeshFunctionGradFE mf_o2_grad_u1(fes_o1_ptr, coeff_vec_u1);
    const lf::fe::MeshFunctionGradFE mf_o2_grad_u2(fes_o1_ptr, coeff_vec_u2);

    // Exact solution for the first component of the velocity
    auto u1 = [&v_ex](Eigen::Vector2d x) -> double { return v_ex(x)[0]; };
    const lf::mesh::utils::MeshFunctionGlobal mf_u1{u1};
    // Exact solution for the gradient of $v_1$
    const lf::mesh::utils::MeshFunctionGlobal mf_grad_u1{grad_v1};
    // Exact solution second component of  the velocity
    auto u2 = [&v_ex](Eigen::Vector2d x) -> double { return v_ex(x)[1]; };
    const lf::mesh::utils::MeshFunctionGlobal mf_u2{u2};
    // Exact solution for the gradient of $v_2$
    const lf::mesh::utils::MeshFunctionGlobal mf_grad_u2{grad_v2};
    // Mesh function for exact solution pressure
    const lf::mesh::utils::MeshFunctionGlobal mf_p{p_ex};
    // compute errors with 5th order quadrature rules
    L2err_u1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_u1 - mf_u1), 4));
    L2err_u2 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_u2 - mf_u2), 4));
    H1err_u1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_grad_u1 - mf_grad_u1),
        4));
    H1err_u2 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o2_grad_u2 - mf_grad_u2),
        4));
    L2err_p = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_o1_p - mf_p), 4));
    errs.emplace_back(dofh.NumDofs(), L2err_u1, L2err_u2, H1err_u1, H1err_u2,
                      L2err_p);
  }
  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
  std::cout.precision(3);
  std::cout << std::endl
            << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 err(v1)" << std::setw(16) << "L2 err(v2)" << std::setw(16)
            << "H1 err(v1)" << std::setw(16) << "H1 err(v2)" << std::setw(16)
            << "L2 err(p)" << '\n';
  std::cout << "---------------------------------------------" << '\n';
  for (const auto& err : errs) {
    auto [N, L2err_u1, L2err_u2, H1err_u1, H1err_u2, L2err_p] = err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << L2err_u1 << std::setw(16) << L2err_u2 << std::setw(16)
             << H1err_u1 << std::setw(16) << H1err_u2 << std::setw(16)
             << L2err_p << '\n';
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << L2err_u1 << std::setw(16) << L2err_u2 << std::setw(16)
              << H1err_u1 << std::setw(16) << H1err_u2 << std::setw(16)
              << L2err_p << '\n';
  }
}

}  // namespace StokesStabP1FEM
