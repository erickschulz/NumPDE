/**
 * @file magstat2d.cc
 * @brief NPDE homework MagStat2D code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#include "magstat2d.h"

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_constant.h>

namespace MagStat2D {

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

void testCvgMagStat2DWhitneyFEM(unsigned int refsteps) {
  using namespace std::numbers;
  // Constant coefficient mu = 1
  lf::mesh::utils::MeshFunctionConstant<double> mf_one(1.0);
  // ********** Part I: Manufactured solution **********
  // Domain: unit square
  // A divergence-free vector field with zero tangential components
  // on the boundary of unit sqaure
  auto u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-std::cos(pi * x[0]) * std::sin(pi * x[1]),
                           std::sin(pi * x[0]) * std::cos(pi * x[1]));
  };
  // The solution component p is zero in this case
  auto p = [](Eigen::Vector2d /*x*/) -> double { return 0.0; };

  // Right-hand side source current; u is an eigenfunction of the vector
  // Laplacian.
  auto f = [&u](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 2 * pi * pi * u(x);
  };
  // Wrap both vectorfields into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  lf::mesh::utils::MeshFunctionGlobal mf_f(f);
  lf::mesh::utils::MeshFunctionGlobal mf_p(p);
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
    Eigen::VectorXd dofvec = solveMagStat2DBVP(dofh, mf_one, mf_f);
    // Build MeshFunction representing the 1-form component of the solution
    const MeshFunctionWF1 mf_sol_u(dofh, dofvec);
    const MeshFunctionWF0 mf_sol_p(dofh, dofvec);
    // Compute L2 norm of the error
    double L2err_u = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_u - mf_u), 2));
    double L2err_p = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol_p - mf_p), 2));
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

}  // namespace MagStat2D
