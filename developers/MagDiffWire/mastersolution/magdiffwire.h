/**
 * @file magdiffwire.h
 * @brief NPDE homework MagDiffWire code
 * @author Ralf Hiptmair
 * @date June 2025
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef MagDiffWire_H_
#define MagDiffWire_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembly_types.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/fix_dof.h>
#include <lf/base/lf_assert.h>
#include <lf/base/types.h>
#include <lf/fe/fe.h>
#include <lf/fe/loc_comp_ellbvp.h>
#include <lf/fe/scalar_fe_space.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

namespace MagDiffWire {

/** @brief Function p.w. constant on two subdomains
 *
 * The distinction between the two subdomains is based on a flag array index
 * by the cells of the mesh.
 */
class MeshFunctionPWConst {
 public:
  MeshFunctionPWConst(const lf::mesh::utils::CodimMeshDataSet<bool>& flags,
                      double val_true, double val_false)
      : flags_(flags), val_true_(val_true), val_false_(val_false) {}
  // Evaluation operator: returns one value for cells flagged 'true', the other
  // value for cells flagged 'false'
  std::vector<double> operator()(const lf::mesh::Entity& cell,
                                 const Eigen::MatrixXd& local) const;

 private:
  const lf::mesh::utils::CodimMeshDataSet<bool>& flags_;
  double val_true_;
  double val_false_;
};

/** @brief Computation of extended M matrix in triplet format =
 *
 */
std::pair<lf::assemble::COOMatrix<double>, lf::assemble::COOMatrix<double>>
buildExtMOLMatrices(std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
                    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags,
                    double sigma_c, double mu_c);

/* @brief SDIRK-2 timestepping for magnetic diffusion equation */
/* SAM_LISTING_BEGIN_3 */
template <typename SCALARFUNCTOR, typename RECORDER = std::function<
                                      void(double, const Eigen::VectorXd&)>>
void sdirkMagDiffWire(
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags, double sigma_c,
    double mu_c, SCALARFUNCTOR I_source, unsigned int M, double T_final,
    Eigen::VectorXd& mu_vec,
    RECORDER rec = [](double, const Eigen::VectorXd&) -> void {}) {
  std::cout << "sdirMagDiffWire_org\n";
  // Obtain extended Galerkin matrices as lf::asseemble::COOMatrix
  auto [Mt, At] = buildExtMOLMatrices(fes_p, Oc_flags, sigma_c, mu_c);
  const size_t N = At.cols() - 1;
  LF_ASSERT_MSG(At.rows() == N + 1, "A must be square");
  LF_ASSERT_MSG(Mt.rows() == N + 1, "M wrong row number");
  LF_ASSERT_MSG(Mt.cols() == N + 1, "M wrong col number");
  LF_ASSERT_MSG(mu_vec.size() == N + 1, "State vector must have N+1 componnts");
  // We have built a sparse matrix containing the diagonal block(s) of the
  // SDIRK2 linear system, because this has to be supplied to Eigen's built-in
  // sparse elimination solver.
  auto& At_triplet_vec{At.triplets()};
  auto& Mt_triplet_vec{Mt.triplets()};
  // Form weighted sum of sparse matrices in triplet format by merging triplet
  // vectors
  const double zeta = 1.0 - 0.5 * std::sqrt(2.0);  // SDIRK parameter
  const double tau = T_final / M;
  std::vector<Eigen::Triplet<double>> Dt_triplet_vec{};
  for (const auto& triplet : Mt_triplet_vec) {
    Dt_triplet_vec.push_back(triplet);
  }
  for (const auto& triplet : At_triplet_vec) {
    Dt_triplet_vec.emplace_back(triplet.row(), triplet.col(),
                                zeta * tau * triplet.value());
  }
  // Convert to CRS format
  Eigen::SparseMatrix<double> D(N + 1, N + 1);
  D.setFromTriplets(Dt_triplet_vec.begin(), Dt_triplet_vec.end());
  D.makeCompressed();
  // For the sake of efficiency: LU-decomposition outside the actual
  // timestepping loop
  Eigen::SparseLU<Eigen::SparseMatrix<double>> DLU;
  DLU.compute(D);
  LF_VERIFY_MSG(DLU.info() == Eigen::Success, "LU decomposition failed");

// Main timestepping loop
// A few auxliary vectors
#if SOLUTION
  // Right hand side vector
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N + 1);
  // Vector holding M * current state
  Eigen::VectorXd Mmu(N + 1);
  // Temporary vector
  Eigen::VectorXd tmp(N + 1);
#else
/* **********************************************************************
   Define auxiliary vectors, if any, here
   ********************************************************************** */
#endif
  // Current time $t_{k-1}$
  double t = 0.0;
  for (int k = 1; k <= M; ++k, t += tau) {
#if SOLUTION
    rhs = Mt.MatVecMult(1.0, mu_vec);
    const double rhsN = rhs[N];
    rhs[N] += zeta * tau * I_source(t + zeta * tau);
    tmp = DLU.solve(rhs);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    rhs[N] = rhsN + tau * I_source(t + tau);
    At.MatVecMult(-tau * (1 - zeta), tmp, rhs);
    mu_vec = DLU.solve(rhs);
#else
    /* **********************************************************************
       Your loop body code here
       ********************************************************************** */
#endif
    rec(t, mu_vec);
  }
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename SCALARFUNCTOR, typename RECORDER = std::function<
                                      void(double, const Eigen::VectorXd&)>>
void sdirkMagDiffWire_alt(
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags, double sigma_c,
    double mu_c, SCALARFUNCTOR I_source, unsigned int M, double T_final,
    Eigen::VectorXd& mu_vec,
    RECORDER rec = [](double, const Eigen::VectorXd&) -> void {}) {
  std::cout << "sdirMagDiffWire_alt\n";
  // Obtain extended Galerkin matrices as lf::asseemble::COOMatrix
  auto [Mt, At] = buildExtMOLMatrices(fes_p, Oc_flags, sigma_c, mu_c);
  const size_t N = At.cols() - 1;
  LF_ASSERT_MSG(At.rows() == N + 1, "A must be square");
  LF_ASSERT_MSG(Mt.rows() == N + 1, "M wrong row number");
  LF_ASSERT_MSG(Mt.cols() == N + 1, "M wrong col number");
  LF_ASSERT_MSG(mu_vec.size() == N + 1, "State vector must have N+1 componnts");
  // We have built a sparse matrix containing the diagonal block(s) of the
  // SDIRK2 linear system, because this has to be supplied to Eigen's built-in
  // sparse elimination solver.
  auto& At_triplet_vec{At.triplets()};
  auto& Mt_triplet_vec{Mt.triplets()};
  // Form weighted sum of sparse matrices in triplet format by merging triplet
  // vectors
  const double zeta = 1.0 - 0.5 * std::sqrt(2.0);  // SDIRK parameter
  const double tau = T_final / M;
  std::vector<Eigen::Triplet<double>> Dt_triplet_vec{};
  for (const auto& triplet : Mt_triplet_vec) {
    Dt_triplet_vec.push_back(triplet);
  }
  for (const auto& triplet : At_triplet_vec) {
    Dt_triplet_vec.emplace_back(triplet.row(), triplet.col(),
                                zeta * tau * triplet.value());
  }
  // Convert to CRS format
  Eigen::SparseMatrix<double> D(N + 1, N + 1);
  D.setFromTriplets(Dt_triplet_vec.begin(), Dt_triplet_vec.end());
  D.makeCompressed();
  // For the sake of efficiency: LU-decomposition outside the actual
  // timestepping loop
  Eigen::SparseLU<Eigen::SparseMatrix<double>> DLU;
  DLU.compute(D);
  LF_VERIFY_MSG(DLU.info() == Eigen::Success, "LU decomposition failed");

  // Main timestepping loop, increment-based implementation
  // Current time $t_{k-1}$
  double t = 0.0;
  for (int k = 1; k <= M; ++k, t += tau) {
#if SOLUTION
    const Eigen::VectorXd Amu = At.MatVecMult(1.0, mu_vec);
    Eigen::VectorXd rho_vec = Eigen::VectorXd::Zero(N + 1);
    rho_vec[N] = I_source(t + zeta * tau);
    const Eigen::VectorXd kappa1 = DLU.solve(rho_vec - Amu);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    rho_vec[N] = I_source(t + tau);
    const Eigen::VectorXd kappa2 = DLU.solve(
        rho_vec - tau * (1 - zeta) * At.MatVecMult(1.0, kappa1) - Amu);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    mu_vec += tau * (1 - zeta) * kappa1 + tau * zeta * kappa2;
#else
    /* **********************************************************************
       Your loop body code here
       ********************************************************************** */
#endif
    rec(t, mu_vec);
  }
}
/* SAM_LISTING_END_4 */

#if SOLUTION

/* SAM_LISTING_BEGIN_5 */
template <typename SCALARFUNCTOR, typename RECORDER = std::function<
                                      void(double, const Eigen::VectorXd&)>>
void sdirkMagDiffWire_mat(
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags, double sigma_c,
    double mu_c, SCALARFUNCTOR I_source, unsigned int M, double T_final,
    Eigen::VectorXd& mu_vec,
    RECORDER rec = [](double, const Eigen::VectorXd&) -> void {}) {
  std::cout << "sdirkMagDiffWire_mat\n";
  // Obtain extended Galerkin matrices as lf::asseemble::COOMatrix
  auto [Mt, At] = buildExtMOLMatrices(fes_p, Oc_flags, sigma_c, mu_c);
  const size_t N = At.cols() - 1;
  LF_ASSERT_MSG(At.rows() == N + 1, "A must be square");
  LF_ASSERT_MSG(Mt.rows() == N + 1, "M wrong row number");
  LF_ASSERT_MSG(Mt.cols() == N + 1, "M wrong col number");
  LF_ASSERT_MSG(mu_vec.size() == N + 1, "State vector must have N+1 componnts");
  // Construct sparse matrices
  const Eigen::SparseMatrix<double> A_mat{At.makeSparse()};
  const Eigen::SparseMatrix<double> M_mat{Mt.makeSparse()};
  // Parameters for SDIRK-2 timestepping
  const double zeta = 1.0 - 0.5 * std::sqrt(2.0);  // SDIRK parameter
  const double tau = T_final / M;
  const Eigen::SparseMatrix<double> D_mat = M_mat + tau * zeta * A_mat;
  // For the sake of efficiency: LU-decomposition outside the actual
  // timestepping loop
  Eigen::SparseLU<Eigen::SparseMatrix<double>> DLU;
  DLU.compute(D_mat);
  LF_VERIFY_MSG(DLU.info() == Eigen::Success, "LU decomposition failed");

  // Main timestepping loop
  // Current time $t_{k-1}$
  double t = 0.0;
  for (int k = 1; k <= M; ++k, t += tau) {
    const Eigen::VectorXd Amu = A_mat * mu_vec;
    Eigen::VectorXd rho_vec = Eigen::VectorXd::Zero(N + 1);
    rho_vec[N] = I_source(t + zeta * tau);
    const Eigen::VectorXd kappa1 = DLU.solve(rho_vec - Amu);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    rho_vec[N] = I_source(t + tau);
    const Eigen::VectorXd kappa2 =
        DLU.solve(rho_vec - tau * (1 - zeta) * A_mat * kappa1 - Amu);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    mu_vec += tau * (1 - zeta) * kappa1 + tau * zeta * kappa2;
    rec(t, mu_vec);
  }
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
template <typename SCALARFUNCTOR, typename RECORDER = std::function<
                                      void(double, const Eigen::VectorXd&)>>
void sdirkMagDiffWire_new(
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fes_p,
    lf::mesh::utils::CodimMeshDataSet<bool> Oc_flags, double sigma_c,
    double mu_c, SCALARFUNCTOR I_source, unsigned int M, double T_final,
    Eigen::VectorXd& mu_vec,
    RECORDER rec = [](double, const Eigen::VectorXd&) -> void {}) {
  std::cout << "sdirkMagDiffWire_new\n";
  // Obtain extended Galerkin matrices as lf::asseemble::COOMatrix
  auto [Mt, At] = buildExtMOLMatrices(fes_p, Oc_flags, sigma_c, mu_c);
  const size_t N = At.cols() - 1;
  LF_ASSERT_MSG(At.rows() == N + 1, "A must be square");
  LF_ASSERT_MSG(Mt.rows() == N + 1, "M wrong row number");
  LF_ASSERT_MSG(Mt.cols() == N + 1, "M wrong col number");
  LF_ASSERT_MSG(mu_vec.size() == N + 1, "State vector must have N+1 componnts");
  // Construct sparse matrices
  const Eigen::SparseMatrix<double> A_mat{At.makeSparse()};
  const Eigen::SparseMatrix<double> M_mat{Mt.makeSparse()};
  // Parameters for SDIRK-2 timestepping
  const double zeta = 1.0 - 0.5 * std::sqrt(2.0);  // SDIRK parameter
  const double tau = T_final / M;
  const Eigen::SparseMatrix<double> D_mat = M_mat + tau * zeta * A_mat;
  // For the sake of efficiency: LU-decomposition outside the actual
  // timestepping loop
  Eigen::SparseLU<Eigen::SparseMatrix<double>> DLU;
  DLU.compute(D_mat);
  LF_VERIFY_MSG(DLU.info() == Eigen::Success, "LU decomposition failed");

  // Main timestepping loop
  // Current time $t_{k-1}$
  double t = 0.0;
  for (int k = 1; k <= M; ++k, t += tau) {
    const Eigen::VectorXd Mmu = M_mat * mu_vec;
    Eigen::VectorXd rho_vec = Eigen::VectorXd::Zero(N + 1);
    rho_vec[N] = I_source(t + zeta * tau);
    const Eigen::VectorXd gamma1 = DLU.solve(Mmu + tau * zeta * rho_vec);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    rho_vec[N] = I_source(t + tau);
    const Eigen::VectorXd gamma2 =
        DLU.solve(Mmu + tau * rho_vec - tau * (1 - zeta) * A_mat * gamma1);
    LF_VERIFY_MSG(DLU.info() == Eigen::Success, "Solving LSE failed");
    mu_vec = gamma2;
    rec(t, mu_vec);
  }
}
/* SAM_LISTING_END_6 */

#endif

/** @brief Convergence test with manufactured solution */
void testCvgMagDiffWire(unsigned int refsteps);

}  // namespace MagDiffWire

#endif
