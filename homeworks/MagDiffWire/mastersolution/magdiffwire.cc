/**
 * @file magdiffwire.cc
 * @brief NPDE homework MagDiffWire code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#include "magdiffwire.h"

#include <lf/fe/fe_tools.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/mesh_function_unary.h>
#include <lf/uscalfe/fe_space_lagrange_o2.h>

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <numbers>

namespace MagDiffWire {

std::vector<double> MeshFunctionPWConst::operator()(
    const lf::mesh::Entity &cell, const Eigen::MatrixXd &local) const {
  const size_t n_pts = local.cols();
  if (flags_(cell)) {
    return std::vector<double>(n_pts, val_true_);
  }
  return std::vector<double>(n_pts, val_false_);
}

/* SAM_LISTING_BEGIN_1 */
std::pair<lf::assemble::COOMatrix<double>, lf::assemble::COOMatrix<double>>
buildExtMOLMatrices(std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
                    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags,
                    double sigma_c, double mu_c) {
  // Obtain handle to dof handler and mesh
  const lf::assemble::DofHandler &dofh{fes_p->LocGlobMap()};
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p = dofh.Mesh();
  // Dimension of finite element space
  const lf::base::size_type N(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> Mt(N + 1, N + 1);
  lf::assemble::COOMatrix<double> At(N + 1, N + 1);
  // Mesh function providing conductivity
  MeshFunctionPWConst mf_sigma(Oc_flags, sigma_c, 0.0);
  // Mesh function providing permeability mu
  LF_ASSERT_MSG(mu_c > 0, "mu must be positive");
  MeshFunctionPWConst mf_mu(Oc_flags, 1.0 / mu_c, 1.0);
  // ELEMENT\_MATRIX\_PROVIDERs: LehrFEM utility classes
  lf::fe::MassElementMatrixProvider emp_M(fes_p, mf_sigma);
  lf::fe::DiffusionElementMatrixProvider emp_A(fes_p, mf_mu);
  // Invoke assembly on cells (co-dimension = 0)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, emp_M, Mt);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, emp_A, At);

  // Fill last row and column of At; exploit partition of unity property of
  // global shape functions: the vector c is the product of M and a vector of
  // all ones.
  Eigen::VectorXd vec_one = Eigen::VectorXd::Constant(N + 1, 1.0);
  vec_one[N] = 0.0;
  const Eigen::VectorXd c = Mt.MatVecMult(1.0, vec_one);
  for (int j = 0; j < N; ++j) {
    At.AddToEntry(N, j, c[j]);
    At.AddToEntry(j, N, c[j]);
  }
  // Suppress basis functions on the boundary
  // Flag \cor{any} entity located on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
  // Flag vector for d.o.f. on the boundary
  std::vector<std::pair<bool, double>> ess_dof_select(N + 1, {false, 0.0});
  // Run through all d.o.f.s and check whether they are associated with a mesh
  // entity on the boundary. Then mark those d.o.f.s in the flag array
  for (lf::assemble::gdof_idx_t j = 0; j < N; ++j) {
    if (bd_flags(dofh.Entity(j))) {
      ess_dof_select[j].first = true;
      ess_dof_select[j].second = 0.0;
    }
  }
  // Here we have to use FIxFlaggedSolutionComponents, because we have to
  // decouple d.o.f.s on the boundary completely from the other d.o.f.s.
  // FixFlaggedSolutionCompAlt would not achieve this
  Eigen::VectorXd dummy{Eigen::VectorXd::Zero(N + 1)};
  // Modify the matrices as in  \lref{eq:ffsc}
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      At, dummy);
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_dof_select](lf::assemble::glb_idx_t dof_idx)
          -> std::pair<bool, double> { return ess_dof_select[dof_idx]; },
      Mt, dummy);
  return {Mt, At};
}
/* SAM_LISTING_END_1 */

void testCvgMagDiffWire(unsigned int refsteps) {
  using namespace std::numbers;
  using std::cos;
  using std::exp;
  using std::sin;
  // Manufactured solution on the unit square
  auto u = [](Eigen::VectorXd x, double t) -> double {
    return sin(pi * x[0]) * sin(pi * x[1]) * exp(-2 * pi * pi * t);
  };
  auto u0 = [&u](Eigen::VectorXd x) -> double { return u(x, 0.0); };
  auto grad_u = [](Eigen::VectorXd x, double t) -> Eigen::Vector2d {
    return pi * exp(-2 * pi * pi * t) *
           Eigen::Vector2d(cos(pi * x[0]) * sin(pi * x[1]),
                           sin(pi * x[0]) * cos(pi * x[1]));
  };
  lf::mesh::utils::MeshFunctionGlobal mf_u0(u0);
  const double T_final = 0.5;  // Final time T
                               // Solution at final time
  auto u_fin = [&u, T_final](Eigen::Vector2d x) -> double {
    return u(x, T_final);
  };
  auto grad_u_fin = [&grad_u, T_final](Eigen::VectorXd x) -> Eigen::Vector2d {
    return grad_u(x, T_final);
  };
  // Wrap into MeshFunctions
  lf::mesh::utils::MeshFunctionGlobal mf_u_fin(u_fin);
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u_fin(grad_u_fin);

  // Exciting current
  auto I_source = [](double t) -> double {
    return 4.0 / (pi * pi) * exp(-2 * pi * pi * t);
  };

  // Generate a small unstructured triangular mesh of the unit square
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Generate sequence of meshes by uniform regular refinement
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy &multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);

  // Number of levels
  const int L = multi_mesh.NumLevels();
  // Loop over the levels
  const unsigned int M0 = 20;  // No. of timesteps on coarsest mesh
  unsigned int M = M0;         // Number of timesteps
  // Table of various error norms
  std::vector<std::tuple<size_t, double, double>> errs;
  for (int level = 0; level < L; ++level, M *= 2) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Initialize Lagrangian finite element space
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO2<double>> fes_o2_ptr =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);
    // Number of FE degrees of freedom
    const size_t N = fes_o2_ptr->LocGlobMap().NumDofs();
    // Setting: constant coefficients = 1, all cells flagged as located in the
    // conducting domain
    const lf::mesh::utils::CodimMeshDataSet<bool> all_cells(lev_mesh_p, 0,
                                                            true);
    // Set initial data: nodal interpolant of u(x,0)
    Eigen::VectorXd mu_vec(N + 1);
    mu_vec.head(N) = lf::fe::NodalProjection(*fes_o2_ptr, mf_u0);
    mu_vec[N] = 0.0;
    // Main timestepping loop
    sdirkMagDiffWire(fes_o2_ptr, all_cells, 1.0, 1.0, I_source, M, T_final,
                     mu_vec);
    // Compute errors at final time
    const lf::fe::MeshFunctionFE mf_uh_fin(fes_o2_ptr, mu_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_uh_fin(fes_o2_ptr, mu_vec);
    const double L2_err = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_uh_fin - mf_u_fin), 4));
    const double H1s_err = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p,
        lf::mesh::utils::squaredNorm(mf_grad_uh_fin - mf_grad_u_fin), 4));
    errs.emplace_back(N, L2_err, H1s_err);
  }
  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
  std::cout.precision(3);
  std::cout << std::endl
            << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 err" << std::setw(16) << "H1s err" << '\n';
  std::cout << "---------------------------------------------" << '\n';
  for (const auto &err : errs) {
    auto [N, L2err, H1err] = err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << L2err << std::setw(16) << H1err << '\n';
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << L2err << std::setw(16) << H1err << '\n';
  }
}

}  // namespace MagDiffWire
