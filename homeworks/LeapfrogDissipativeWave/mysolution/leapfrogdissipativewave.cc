/**
 * @file leapfrogdissipativewave.cc
 * @brief NPDE homework LeapfrogDissipativeWave code
 * @author R. Hiptmair
 * @date July 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include "leapfrogdissipativewave.h"

#include <math.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include "lf/mesh/test_utils/test_meshes.h"

namespace LeapfrogDissipativeWave {

/* SAM_LISTING_BEGIN_1 */
std::tuple<lf::assemble::COOMatrix<double>, lf::assemble::COOMatrix<double>,
           lf::assemble::COOMatrix<double>>
computeGalerkinMatrices(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space_p) {
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Instantiating N_dofs x N_dofs zero matrices in triplet format
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> B_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  /* Creating dummy coefficient-functions as Lehrfem++ mesh functions */
  // These coefficient-functions are needed in the class templates
  // ReactionDiffusionElementMatrixProvider and MassEdgeMatrixProvider
  auto zero_mf = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto one_mf = lf::mesh::utils::MeshFunctionConstant(1.0);
  /* Initialization of local matrices builders */
  // Initialize objects taking care of local computations for volume integrals
  lf::uscalfe::ReactionDiffusionElementMatrixProvider M_locmat_builder(
      fe_space_p, zero_mf, one_mf);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider A_locmat_builder(
      fe_space_p, one_mf, zero_mf);
  // Initialize objects taking care of local computations for boundary integrals
  // Creating a predicate that will guarantee that the computations for the
  // boundary mass matrix are carried only on the edges of the mesh
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  lf::uscalfe::MassEdgeMatrixProvider B_locmat_builder(fe_space_p, one_mf,
                                                       bd_flags);
  /* Assembling the Galerkin matrices */
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix.
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_locmat_builder, M_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A_locmat_builder, A_COO);
  // Invoke assembly on edges (co-dimension = 1 as first argument)
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, B_locmat_builder, B_COO);
  return {M_COO, B_COO, A_COO};
}
/* SAM_LISTING_END_1 */

double timestepScalarTestProblem(unsigned int M) {
  // Initial values
  const double u0 = 1.0;
  const double v0 = -0.5;
  // Tiemstep size
  const double tau = 1.0 / M;
  // Approximations of solution and derivatives
  double u = u0;
  // Special initial step
  double v = v0 - 0.5 * tau * (u0 + v0);
  // Timstepping loop
  for (int j = 1; j <= M; ++j) {
    u = u + tau * v;
    v = ((1.0 - 0.5 * tau) * v - tau * u) / (1.0 + 0.5 * tau);
  }
  return u;
}

void convergenceScalarTestProblem(void) {
  std::cout << "Convergence test for scalar model problem" << std::endl;
  // Exact solution of y''+y'+y = 0
  auto u_exact = [](double t) -> double {
    return std::exp(-t / 2.0) * std::cos(0.5 * std::sqrt(3) * t);
  };
  double uT = u_exact(1.0);

  const int nIter = 13;  // total number of runs
  unsigned int M;        // number of equidistant steps
  double tau;            // time step "tau"
  // Error between the approx solutions as given by the implicit method
  // and the exact solution vector computed from the anlytic formula vector
  std::vector<double> errs(nIter);
  unsigned int M_stored[13] = {10,  20,  30,  40,  50,  60, 80,
                               100, 160, 200, 320, 500, 640};
  for (int k = 0; k < nIter; k++) {
    // Number of timesteps
    M = M_stored[k];
    errs[k] = std::abs(uT - timestepScalarTestProblem(M));
  }
  std::cout << std::left << std::fixed << std::setprecision(16) << std::setw(10)
            << "M" << std::setw(20) << "error" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << std::setw(10) << M_stored[k] << std::setw(20) << errs[k]
              << std::endl;
  }
}

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd timestepDissipativeWaveEquation(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space_p,
    double T, unsigned int M, Eigen::VectorXd mu0, Eigen::VectorXd nu0) {
  // Pointer and reference to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Obtain local->global index mapping
  // for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space = number of nodes of the mesh
  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_VERIFY_MSG(mu0.size() == N_dofs, "Wrong length of mu0");
  LF_VERIFY_MSG(nu0.size() == N_dofs, "Wrong length of nu0");
  // The solution vector at time T
  Eigen::VectorXd mu = mu0;

  // ========================================
  // Your code here
  // ========================================
  return mu;
}
/* SAM_LISTING_END_2 */

void convergenceDissipativeLeapfrog(unsigned int reflevels, double T,
                                    unsigned int M0, unsigned int Mfac) {
  // Setting for convergence test: Initial data as mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_u0{[](Eigen::Vector2d x) -> double {
    return std::cos(M_PI * (2 * x[0] - 1)) * std::cos(M_PI * (2 * x[1] - 1));
  }};
  lf::mesh::utils::MeshFunctionConstant mf_v0{1.0};

  // Generate a coarse mesh of a unit-square domain
  const int selector = 0;
  const double scale = 1.0 / 3.0;
  std::shared_ptr<lf::mesh::Mesh> c_mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector, scale);
  // Obtain a pointer to a hierarchy of nested meshes
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(c_mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput information about hierarchy of nested meshes
  std::cout << "\t Unit square: Sequence of nested meshes used for convergence "
               "test\n";
  multi_mesh.PrintInfo(std::cout);
  // Number of levels
  lf::base::size_type L = multi_mesh.NumLevels();
  // Instantiate FE spaces for all levels and store pointers to them
  std::vector<std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>>>
      fe_space_ptrs;
  for (lf::base::size_type level = 0; level < L; ++level) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);
    fe_space_ptrs.push_back(
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p));
  }
  // LEVEL LOOP: Do computations on all levels
  // Solution coefficient vectors
  std::vector<Eigen::VectorXd> mu_vecs;
  // Vector for storing norms of differences
  std::vector<std::tuple<lf::base::size_type, lf::base::size_type, double,
                         double, double, double>>
      norms{};
  unsigned int M = M0;
  for (lf::base::size_type level = 0; level < L; ++level, M *= Mfac) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = multi_mesh.getMesh(level);
    const lf::mesh::Mesh& mesh{*mesh_p};
    // Set up global FE space; second-order Lagrangian finite elements
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space_p =
        fe_space_ptrs[level];
    // Interpolate initial data into FE space
    const Eigen::VectorXd u0_vec = lf::fe::NodalProjection(*fe_space_p, mf_u0);
    const Eigen::VectorXd v0_vec = lf::fe::NodalProjection(*fe_space_p, mf_v0);
    // Coefficient vector of solution at final time
    mu_vecs.push_back(
        timestepDissipativeWaveEquation(fe_space_p, T, M, u0_vec, v0_vec));
  }
  const lf::mesh::Mesh& fine_mesh = *multi_mesh.getMesh(L - 1);
  // create mesh functions representing solution / gradient of solution
  const lf::fe::MeshFunctionFE mf_sol(fe_space_ptrs[L - 1], mu_vecs[L - 1]);
  const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space_ptrs[L - 1],
                                               mu_vecs[L - 1]);
  M = M0;
  for (lf::base::size_type level = 0; level < L - 1; ++level, M *= Mfac) {
    // Transfer the coarse solution onto the finest mesh
    const lf::fe::MeshFunctionFE mf_coarse(fe_space_ptrs[level],
                                           mu_vecs[level]);
    const lf::refinement::MeshFunctionTransfer mf_coarse_sol(
        multi_mesh, mf_coarse, level, L - 1);
    const lf::fe::MeshFunctionGradFE mf_grad_coarse(fe_space_ptrs[level],
                                                    mu_vecs[level]);
    const lf::refinement::MeshFunctionTransfer mf_grad_coarse_sol(
        multi_mesh, mf_grad_coarse, level, L - 1);
    // Compute norms
    const double L2norm = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        fine_mesh, lf::mesh::utils::squaredNorm(mf_coarse_sol), 4));
    const double H1snorm = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        fine_mesh, lf::mesh::utils::squaredNorm(mf_grad_coarse_sol), 4));
    // Norms of differences of solutions on current level and previous level
    const double L2n_diff = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        fine_mesh, lf::mesh::utils::squaredNorm(mf_sol - mf_coarse_sol), 4));
    const double H1sn_diff = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        fine_mesh,
        lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_coarse_sol), 4));
    // Number of degrees of freedom
    const lf::base::size_type N_dofs(
        (fe_space_ptrs[level]->LocGlobMap()).NumDofs());
    std::cout << "level = " << level << ", N_dofs = " << N_dofs << ", M = " << M
              << ", L2norm = " << L2norm << ",  H1snorm = " << H1snorm
              << ", L2n_diff = " << L2n_diff << ",  H1sn_diff = " << H1sn_diff
              << std::endl;
    norms.emplace_back(level, N_dofs, L2norm, H1snorm, L2n_diff, H1sn_diff);
  }
}

void testDissipativeLeapfrog(void) {
  std::cout << "Convergence of timestepping" << std::endl;
  // Generate a triangular test mesh of the unit square
  const int selector = 3;  // Mesh comprising only triangles
  const double scale = 1.0 / 3.0;
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector, scale);
  const lf::mesh::Mesh& mesh{*mesh_p};
  // Set up global FE space; quadratic Lagrangian finite elements
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  // Number of degrees of freedom
  const lf::base::size_type N_dofs((fe_space_p->LocGlobMap()).NumDofs());
  // Initial vectors
  Eigen::VectorXd mu0 = Eigen::VectorXd::LinSpaced(N_dofs, 0.0, 1.0);
  Eigen::VectorXd nu0 = Eigen::VectorXd::Zero(N_dofs);

  // Compute solution at final time for different numbers of timesteps
  // increasing in geometric progression
  unsigned int M = 10;
  unsigned int L = 8;
  std::vector<Eigen::VectorXd> mu_vecs;
  for (int l = 0; l < L; ++l, M *= 2) {
    mu_vecs.push_back(
        timestepDissipativeWaveEquation(fe_space_p, 1.0, M, mu0, nu0));
    std::cout << "M = " << M << ": |mu| = " << mu_vecs.back().norm()
              << std::endl;
  }

  std::cout << "Convergence (L2 norm and H1 semi-norm) w.r.t. to M = : "
            << M / 2 << std::endl;
  M = 10;
  for (int l = 0; l < L - 1; ++l, M *= 2) {
    const lf::fe::MeshFunctionFE mf_sol(fe_space_p, mu_vecs.back());
    const lf::fe::MeshFunctionFE mf_coarse_sol(fe_space_p, mu_vecs[l]);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space_p, mu_vecs.back());
    const lf::fe::MeshFunctionGradFE mf_grad_coarse_sol(fe_space_p, mu_vecs[l]);
    const double L2n = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        mesh, lf::mesh::utils::squaredNorm(mf_sol - mf_coarse_sol), 4));
    const double H1sn_diff = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        mesh, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_coarse_sol),
        4));
    std::cout << " M = " << M << "; diff L2 norm = " << L2n
              << "; diff H1 semi-norm = " << H1sn_diff << std::endl;
  }
}

}  // namespace LeapfrogDissipativeWave
