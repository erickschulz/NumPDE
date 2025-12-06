/**
 * @file hodgelaplacian2d.cc
 * @brief NPDE homework HodgeLaplacian2D code
 * @author Ralf Hoptmair
 * @date May 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#include "hodgelaplacian2d.h"

#include <lf/assemble/assembly_types.h>
#include <lf/base/lf_assert.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_global.h>

#include <cstddef>
#include <iomanip>
#include <memory>
#include <span>

namespace HodgeLaplacian2D {
/* SAM_LISTING_BEGIN_1 */
HodgeLaplacian2DElementMatrixProvider::ElemMat
HodgeLaplacian2DElementMatrixProvider::Eval_ref(const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  LF_VERIFY_MSG(cell.Geometry()->isAffine(),
                "Triangle must have straight edges");
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Compute gradients of barycentric coordinate functions, see
  // \lref{cpp:gradbarycoords}. Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns
  const auto G{X.inverse().block<2, 3>(1, 0)};
#if SOLUTION
  // Dummy lambda functions for barycentric coordinates
  const std::array<std::function<double(Eigen::Vector3d)>, 3> lambda{
      [](Eigen::Vector3d c) -> double { return c[0]; },
      [](Eigen::Vector3d c) -> double { return c[1]; },
      [](Eigen::Vector3d c) -> double { return c[2]; }};
  // Lambda functions for local Whitney 1-forms, see \prbcref{eq:lsf1}.
  // Note that $\grad\lambda_i$ is accessed as G.col(i-1).
  const std::array<std::function<Eigen::Vector2d(Eigen::Vector3d)>, 3> beta{
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return c[0] * G.col(1) - c[1] * G.col(0);
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return c[1] * G.col(2) - c[2] * G.col(1);
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return c[2] * G.col(0) - c[0] * G.col(2);
      }};
  // Barycentric coordinates of the midpoints of the edges for
  // use with the third order 3-point edge midpoint quadrature rule
  // \lref{eq:qr2}
  const std::array<Eigen::Vector3d, 3> mp = {Eigen::Vector3d({0.5, 0.5, 0}),
                                             Eigen::Vector3d({0, 0.5, 0.5}),
                                             Eigen::Vector3d({0.5, 0, 0.5})};
  // Compute the element matrix  for $-\Delta$ and $\cob{\Cs^0_1}$.
  // See also \lref{mc:ElementMatrixLaplLFE}.
  // Eigen::Matrix<double, 3, 3> L = area * G.transpose() * G;
  // Fill the four 3x3 blocks of $\VM_K$
  const double A_ent = 3.0 / (area * area);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // Left upper block $-\VC_K$:
      MK_(i, j) = -(lambda[i](mp[0]) * lambda[j](mp[0]) +
                    lambda[i](mp[1]) * lambda[j](mp[1]) +
                    lambda[i](mp[2]) * lambda[j](mp[2]));
      // Upper right block $\VB_K$ and lower left block $\VB_K^{\top}$:
      const double val =
          ((beta[j](mp[0])) + (beta[j](mp[1])) + (beta[j](mp[2])))
              .dot(G.col(i));
      MK_(i, j + 3) = val;
      MK_(j + 3, i) = val;
      // Lower right block
      MK_(i + 3, j + 3) = A_ent;
    }
  }
  // Correct for orientation mismatch, cf. \prbcref{sp:H}.
  auto relor = cell.RelativeOrientations();
  LF_ASSERT_MSG(relor.size() == 3, "Triangle should have 3 edges!?");
  for (int k = 0; k < 3; ++k) {
    if (relor[k] == lf::mesh::Orientation::negative) {
      // Flip sign of 3+k-th rown and column
      MK_.col(3 + k) *= -1.0;
      MK_.row(3 + k) *= -1.0;
    }
  }
  // Finally multiply with with the quadrature weight.
  MK_ *= area / 3.0;
#else
/* **********************************************************************
   Your code here
   ********************************************************************** */
#endif
  return MK_;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
HodgeLaplacian2DElementMatrixProvider::ElemMat
HodgeLaplacian2DElementMatrixProvider::Eval(const lf::mesh::Entity& cell) {
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
#if SOLUTION
  // Compute the element matrix  for $-\Delta$ and $\cob{\Cs^0_1}$.
  // See also \lref{mc:ElementMatrixLaplLFE}.
  Eigen::Matrix<double, 3, 3> L = area * G.transpose() * G;

  // Fill left upper block $\VC_K$, \prbcref{eq:Mk}
  const Eigen::Matrix3d C =
      (Eigen::Matrix3d() << 2, 1, 1, 1, 2, 1, 1, 1, 2).finished();
  MK_.block(0, 0, 3, 3) = -(area / 12.0) * C;
  // Initialize right upper block $\VB_K$ and left lower block $\VB^\top$,
  // see \prbcref{eq:BKdef}.
  Eigen::Matrix3d B;
  B.col(0) = (L.col(1) - L.col(0)) / 3.0;
  B.col(1) = (L.col(2) - L.col(1)) / 3.0;
  B.col(2) = (L.col(0) - L.col(2)) / 3.0;
  MK_.block(0, 3, 3, 3) = B;
  MK_.block(3, 0, 3, 3) = B.transpose();
  // Set lower right block $\VA_K$, see \prbcref{eq:Akrefc}
  MK_.block(3, 3, 3, 3) = Eigen::Matrix3d::Constant(1.0 / area);
  // Correct for orientation mismatch, cf. \prbcref{sp:H}.
  auto relor = cell.RelativeOrientations();
  LF_ASSERT_MSG(relor.size() == 3, "Triangle should have 3 edges!?");
  for (int k = 0; k < 3; ++k) {
    if (relor[k] == lf::mesh::Orientation::negative) {
      // Flip sign of 3+k-th rown and column
      MK_.col(3 + k) *= -1.0;
      MK_.row(3 + k) *= -1.0;
    }
  }
#else
/* **********************************************************************
   Your code here
   ********************************************************************** */
#endif
  return MK_;
}
/* SAM_LISTING_END_2 */

lf::assemble::COOMatrix<double> buildHodgeLaplacianGalerkinMatrix(
    const lf::assemble::DofHandler& dofh) {
  // Total number of FE d.o.f.s
  lf::assemble::size_type N = dofh.NumDofs();
  // Full Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A(N, N);
#if SOLUTION
  // Set up computation of element matrix
  HodgeLaplacian2DElementMatrixProvider hlemp{};
  // Assemble \cor{full} Galerkin matrix for Whitney FEM for
  // Hodge-Laplacian mixed variational problem
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, hlemp, A);
#else
  /* **********************************************************************
     Your code here
     ********************************************************************** */
#endif
  return A;
}

std::vector<Eigen::Vector2d> MeshFunctionWF1::operator()(
    const lf::mesh::Entity& cell, const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported entity type " << cell.RefEl());
  LF_VERIFY_MSG(cell.Geometry()->isAffine(),
                "Triangle must have straight edges");
  /*
  // Compute gradients of barycentric coordinate functions, see
  // \lref{cpp:gradbarycoords}. Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns.
  // Note that $\grad\lambda_i$ is accessed as G.col(i-1).
  const auto G{X.inverse().block<2, 3>(1, 0)};*/
  // Compute gradients of barycentric coordinate functions and store them in the
  // columns of the $2\times 3$-matrix G
  // clang-format off
  const Eigen::MatrixXd dpt = (Eigen::MatrixXd(2, 1) << 0.0, 0.0).finished();
  const Eigen::Matrix<double, 2, 3> G =
      cell.Geometry()->JacobianInverseGramian(dpt).block(0, 0, 2, 2) *
    (Eigen::Matrix<double, 2, 3>(2,3) << -1, 1, 0,
                                         -1, 0, 1).finished();
  // clang-format on

  // Lambda functions in terms of reference coordinates for local Whitney
  // 1-forms, see \prbcref{eq:lsf1}.
  const std::array<std::function<Eigen::Vector2d(Eigen::Vector2d)>, 3> beta{
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return (1.0 - c[0] - c[1]) * G.col(1) - c[0] * G.col(0);
      },
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return c[0] * G.col(2) - c[1] * G.col(1);
      },
      [&G](Eigen::Vector2d c) -> Eigen::Vector2d {
        return c[1] * G.col(0) - (1.0 - c[0] - c[1]) * G.col(2);
      }};
  //  Obtain local d.o.f.s
  std::span<const lf::assemble::gdof_idx_t> locdofs{
      dofh_.GlobalDofIndices(cell)};
  LF_ASSERT_MSG(locdofs.size() == 6, "Six dofs must belong to every cell");
  // The local d.o.f.s with number 3,4,5 are the coefficients for the local
  // Whitney 1-forms (edge basis functions)
  Eigen::Vector3d wf1ldofs(coeffs_[locdofs[3]], coeffs_[locdofs[4]],
                           coeffs_[locdofs[5]]);
  // Correct for orientation mismatch
  auto relor = cell.RelativeOrientations();
  LF_ASSERT_MSG(relor.size() == 3, "Triangle should have 3 edges!?");
  for (int k = 0; k < 3; ++k) {
    if (relor[k] == lf::mesh::Orientation::negative) {
      // Flip sign of coefficient
      wf1ldofs[k] *= -1;
    }
  }
  std::vector<Eigen::Vector2d> res;
  // Run through evaluation points (columns of argument matrix)
  const size_t no_pts = local.cols();
  for (int j = 0; j < no_pts; ++j) {
    res.emplace_back(wf1ldofs[0] * beta[0](local.col(j)) +
                     wf1ldofs[1] * beta[1](local.col(j)) +
                     wf1ldofs[2] * beta[2](local.col(j)));
  }
  return res;
}

std::vector<double> MeshFunctionWF0::operator()(
    const lf::mesh::Entity& cell, const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported entity type " << cell.RefEl());
  //  Obtain local d.o.f.s
  std::span<const lf::assemble::gdof_idx_t> locdofs{
      dofh_.GlobalDofIndices(cell)};
  LF_ASSERT_MSG(locdofs.size() == 6, "Six dofs must belong to every cell");
  // The local d.o.f.s with number 0,1,2 are the coefficients for the local
  // Whtiney 0-forms (barycentric coordinate functions)
  const Eigen::Vector3d wf0dofs(coeffs_[locdofs[0]], coeffs_[locdofs[1]],
                                coeffs_[locdofs[2]]);
  // Run through evaluation points (columns of argument matrix)
  const size_t no_pts = local.cols();
  std::vector<double> res(no_pts);
  for (int j = 0; j < no_pts; ++j) {
    const Eigen::Vector2d xh{local.col(j)};
    res[j] = (wf0dofs[0] * (1.0 - xh[0] - xh[1]) + wf0dofs[1] * xh[0] +
              wf0dofs[2] * xh[1]);
  }
  return res;
}

std::pair<lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d>>
reconstructNodalFields(const lf::assemble::DofHandler& dofh,
                       const Eigen::VectorXd& coeffs) {
  // Obtain mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  const lf::mesh::Mesh& mesh = *mesh_p;
  LF_ASSERT_MSG(dofh.NumDofs() == (mesh.NumEntities(2) + mesh.NumEntities(1)),
                "DofH must manage 1 dof/node and 1 dof/edge");
  // **************************************************
  // Average values of 1-form component in the nodes of the mesh
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> u_vals(
      mesh_p, 2, Eigen::Vector2d::Constant(0.0));
  // Auxliary array counting the numbers of triangles adjacent to a node
  lf::mesh::utils::CodimMeshDataSet<size_t> nb_cnt(mesh_p, 2, 0);
  // Auxiliary MeshFunction representing 1-form in Whitney finite element space
  const MeshFunctionWF1 mf_u(dofh, coeffs);
  // Reference coordinates of vertices of the triangle
  const Eigen::MatrixXd vert_coord =
      (Eigen::MatrixXd(2, 3) << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0).finished();
  // Run through the cells of the mesh
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    LF_VERIFY_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Unsupported entity type " << cell->RefEl());
    // Obtain values of discrete 1-form in the vertices of the mesh
    std::vector<Eigen::Vector2d> wf1_vertex_vals = mf_u(*cell, vert_coord);
    LF_ASSERT_MSG(wf1_vertex_vals.size() == 3, "Need 3 values for 3 vertices");
    // Fetch nodes of the cell
    std::span<const lf::mesh::Entity* const> nodes{cell->SubEntities(2)};
    LF_ASSERT_MSG(nodes.size() == 3, "Triangle must have 3 nodes");
    for (int k = 0; k < nodes.size(); ++k) {
      nb_cnt(*nodes[k])++;
      u_vals(*nodes[k]) += wf1_vertex_vals[k];
    }
  }
  // **************************************************
  // Collect nodal values of 0-form components and compute averages for nodal
  // vectors of 1-form component
  lf::mesh::utils::CodimMeshDataSet<double> p_vals(mesh_p, 2, 0.0);
  // Run through the nodes of the mesh
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    std::span<const lf::assemble::gdof_idx_t> dof_idx{
        dofh.InteriorGlobalDofIndices(*node)};
    LF_ASSERT_MSG(dof_idx.size() == 1, "One d.o.f. per node1");
    p_vals(*node) = coeffs[dof_idx[0]];
    u_vals(*node) /= nb_cnt(*node);
  }

  return {p_vals, u_vals};
}

/** @brief test of convergence based on manufactured solution
 *
 * Computes L2 norm of error of the FEM solution for the 1-form component on a
 * sequence of meshes generated by regular refinement.
 */
void testCvgHLWhitneyFEM(unsigned int refsteps) {
  using namespace std::numbers;
  // ********** Part I: Manufactured solution **********
  // Domain: unit square
  // A curl-free vectorfield with zero normal components at the boundary of the
  // unit square: satisfies all natural b.c. for HL BVP
  auto u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(std::sin(pi * x[0]) * std::cos(pi * x[1]),
                           std::cos(pi * x[0]) * std::sin(pi * x[1]));
  };
  // Right-hand side source vectorfield; u is an eigenfunction of the vector
  // Laplacian.
  auto f = [&u](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 2 * pi * pi * u(x);
  };
  // Divergence of u
  auto divu = [](Eigen::Vector2d x) -> double {
    return 2 * pi * std::cos(pi * x[0]) * std::cos(pi * x[1]);
  };
  // Wrap both vectorfields into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  lf::mesh::utils::MeshFunctionGlobal mf_divu(divu);
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
  // Table of error norms
  std::vector<std::tuple<size_t, double, double>> L2errs;
  for (int level = 0; level < L; ++level) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // ********** Part III: Solving on a single mesh **********
    // Initialize dof handler for monolithic Whitney FEM
    lf::assemble::UniformFEDofHandler dofh(lev_mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 0}});
    const size_t N = dofh.NumDofs();
    LF_ASSERT_MSG(
        (N == lev_mesh_p->NumEntities(2) + lev_mesh_p->NumEntities(1)),
        "No dof mismatch");
    Eigen::VectorXd dofvec = solveHodgeLaplaceBVP(dofh, mf_f);
    // Build MeshFunction representing the 1-form component of the solution
    const HodgeLaplacian2D::MeshFunctionWF1 mf_sol_u(dofh, dofvec);
    const HodgeLaplacian2D::MeshFunctionWF0 mf_sol_p(dofh, dofvec);
    // Compute L2 norm of the error
    double L2err_u = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_u - mf_u), 2));
    double L2err_p = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_p + mf_divu), 2));
    L2errs.push_back({N, L2err_u, L2err_p});
    std::cout << "L2 errors on level " << level << " : of u  = " << L2err_u
              << ", of p = " << L2err_p << std::endl;
  }
  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
  std::cout.precision(3);
  std::cout << std::endl
            << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 err(u)" << std::setw(16) << "L2 err(p)" << '\n';
  std::cout << "---------------------------------------------" << '\n';
  for (const auto& err : L2errs) {
    auto [N, L2err_u, L2err_p] = err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << L2err_u << std::setw(16) << L2err_p << '\n';
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << L2err_u << std::setw(16) << L2err_p << '\n';
  }
}

}  // namespace HodgeLaplacian2D
