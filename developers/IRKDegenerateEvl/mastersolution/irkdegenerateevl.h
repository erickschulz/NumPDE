/**
 * @file irkdegenerateevl.h
 * @brief NPDE homework IRKDegenerateEvl code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef IRKDegenerateEvl_H_
#define IRKDegenerateEvl_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/fe/scalar_fe_space.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>

namespace IRKDegenerateEvl {

/** @brief Assembly of boundary mass matrix in triplet formst
 *
 * @param fes_p pointer to the scalar-valued H^1-conforming finite element space
 * used a trial and test space
 */
[[nodiscard]] lf ::assemble::COOMatrix<double> buildM(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fes_p);

/** @brief Assembly of Galerkin matrix for the negative Laplacian
 *
 * @param fes_p pointer to the scalar-valued H^1-conforming finite element space
 * used a trial and test space
 */
[[nodiscard]] lf ::assemble::COOMatrix<double> buildA(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fes_p);

#define REPLACE 1
#define REPLACE_MAT Eigen::MatrixXd()

/** @brief Two step Radau timestepping for the semidiscrete evolution problem
 *
 * @param M matrix in front of derivative in MOL ODE
 * @param A matrix from MOL ODE
 * @param Ark Butcher matrix for L-stable implicit RK-SSM
 * @param mu0 basis expansion coefficient vector for intial data
 * @param no_ts number of uniform timesteps
 * @param T final time
 * @param rec object for recording sequence of coefficient vector
 */
/* SAM_LISTING_BEGIN_3 */
template <typename MATRIX,
          typename RECORDER = std::function<void(const Eigen::VectorXd&)>>
Eigen::VectorXd timesteppingIRKMOLODE(
    const lf::assemble::COOMatrix<double>& M,
    const lf::assemble::COOMatrix<double>& A, const MATRIX& Ark,
    const Eigen::VectorXd& mu0, unsigned int no_ts, double T,
    RECORDER&& rec = [](const Eigen::VectorXd& /*mu_vec*/) -> void {}) {
  const Eigen::Index s = Ark.cols();
  const Eigen::Index N = A.cols();
  LF_ASSERT_MSG(s == Ark.rows(), "Butcher matrix must be square");
  LF_ASSERT_MSG(N == A.rows(), "Matrix A must be square");
  LF_ASSERT_MSG((N == M.rows()) and (N == M.cols()),
                "Matrix M same size as A!");
  LF_ASSERT_MSG(mu0.size() == N, "Illegal size of mu0!");
  const double tau = T / no_ts;
#if SOLUTION
  // Initialization of system matrix for comnined increment equations
  lf::assemble::COOMatrix<double> ISM(N * s, N * s);
  const lf::assemble::COOMatrix<double>::TripletVec Mtv{M.triplets()};
  const lf::assemble::COOMatrix<double>::TripletVec Atv{A.triplets()};
  for (int i = 0; i < s; ++i) {
    for (int j = 0; j < s; ++j) {
      if (i == j) {
        for (const Eigen::Triplet<double>& M_trp : Mtv) {
          ISM.AddToEntry(M_trp.row() + i * N, M_trp.col() + j * N,
                         M_trp.value());
        }
      }
      for (const Eigen::Triplet<double>& A_trp : Atv) {
        ISM.AddToEntry(A_trp.row() + i * N, A_trp.col() + j * N,
                       tau * Ark(i, j) * A_trp.value());
      }
    }
  }
  // Advance LU decomposition of system matrix for increment equations
  Eigen::SparseMatrix<double> ISM_crs{ISM.makeSparse()};
  Eigen::SparseLU<Eigen::SparseMatrix<double>> ISLU;
  ISLU.compute(ISM_crs);
  LF_VERIFY_MSG(ISLU.info() == Eigen::Success, "LU decomposition failed");
  // Main timestepping loop
  Eigen::VectorXd mu{mu0};
  rec(mu);
  Eigen::VectorXd incs(N * s);
  Eigen::VectorXd Amu(N);
  for (int k = 0; k < no_ts; ++k) {
    // Solve increment system
    Amu = -A.MatVecMult(1.0, mu);
    incs = ISLU.solve(Amu.replicate(s, 1));
    LF_VERIFY_MSG(ISLU.info() == Eigen::Success,
                  "Solving increment system failed");
    mu += tau * Eigen::Map<Eigen::MatrixXd>(incs.data(), N, s) *
          Ark.row(s - 1).transpose();
    rec(mu);
  }
  return mu;
#else
  /* SAM_LISTING_BEGIN_X */
  Eigen::VectorXd REPLACE_VEC(A.cols());

  // Initialization of system matrix for combined increment equations
  lf::assemble::COOMatrix<double> ISM(REPLACE, REPLACE);
  const lf::assemble::COOMatrix<double>::TripletVec Mtv{M.triplets()};
  const lf::assemble::COOMatrix<double>::TripletVec Atv{A.triplets()};
  for (int i = 0; i < REPLACE; ++i) {
    for (int j = 0; j < REPLACE; ++j) {
      if (i == j) {
        for (const Eigen::Triplet<double>& M_trp : Mtv) {
          ISM.AddToEntry(REPLACE, REPLACE, REPLACE);
        }
      }
      for (const Eigen::Triplet<double>& A_trp : Atv) {
        ISM.AddToEntry(REPLACE, REPLACE, REPLACE);
      }
    }
  }

  // Advance LU decomposition of system matrix for increment equations
  Eigen::SparseMatrix<double> ISM_crs{ISM.makeSparse()};
  Eigen::SparseLU<Eigen::SparseMatrix<double>> ISLU;
  ISLU.compute(ISM_crs);

  // Main timestepping loop
  Eigen::VectorXd mu{mu0};
  rec(mu);
  Eigen::VectorXd incs(REPLACE);
  Eigen::VectorXd Amu(N);
  for (int k = 0; k < no_ts; ++k) {
    // clang-format off
    // Solve increment system
    // *********** Uncomment the following line and replace REPLACE with your code
    // Amu = -A.MatVecMult(REPLACE, REPLACE);
    // *********** Uncomment the following line and replace REPLACE with your code
    // incs = ISLU.solve(Amu.replicate(REPLACE, REPLACE));
    LF_VERIFY_MSG(ISLU.info() == Eigen::Success, "Solving increment system failed");
    // *********** Uncomment the following line and replace REPLACE with your code
    // mu += REPLACE * Eigen::Map<Eigen::MatrixXd>(incs.data(), REPLACE, REPLACE) *  REPLACE;
    // clang-format on
    rec(mu);
  }
  return mu;
#endif
}
/* SAM_LISTING_END_3 */

/** @brief Tabulate norms of approximate solution vs timestep *
 *
 * @tparam U0FUNCTOR functor std::function<double(Eigen::Vector2d)
 * @param fes_p pointer to underlying Lagrangian finite element space
 * @param u0 function providing initial data
 * @param no_ts number of timesteps
 * @param T final time
 * @param o output stream (default std::cout)
 *
 * Timestepping based on 2-stage Radau method of oder 3
 */
/* SAM_LISTING_BEGIN_6 */
template <typename U0FUNCTOR>
std::vector<std::pair<double, double>> tabulateSolNorms(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fes_p,
    U0FUNCTOR&& u0, unsigned int no_ts, double T, std::ostream& o = std::cout) {
  // Build MOL ODE matrices
  lf::assemble::COOMatrix<double> M{IRKDegenerateEvl::buildM(fes_p)};
  lf::assemble::COOMatrix<double> A{IRKDegenerateEvl::buildA(fes_p)};
  // Butcher matrix for 2-stage Radau method
  const Eigen::Matrix<double, 2, 2> Ark =
      (Eigen::Matrix<double, 2, 2>() << 5.0 / 12.0, -1.0 / 12.0, 3.0 / 4.0,
       1.0 / 4.0)
          .finished();
  // Coefficient vector of for initial data
  const Eigen::VectorXd mu0 =
      lf::fe::NodalProjection(*fes_p, lf::mesh::utils::MeshFunctionGlobal(u0));
  // Vector for storing the two desired norms
  // First: L2-norm on the boundary, Second: H1 seminorm on the domain
  std::vector<std::pair<double, double>> it_norms{};
  // "Recorder" lambda function
#if SOLUTION
  auto rec = [&A, &M, &it_norms](const Eigen::VectorXd& mu_vec) -> void {
    it_norms.emplace_back(std::sqrt(mu_vec.dot(M.MatVecMult(1.0, mu_vec))),
                          std::sqrt(mu_vec.dot(A.MatVecMult(1.0, mu_vec))));
  };
#else
  /* **************************************************
   * Your code here
   * The next line is just a dummy implementation
   * ************************************************* */
  auto rec = [](const Eigen::VectorXd&) -> void {};
#endif
  (void)timesteppingIRKMOLODE(M, A, Ark, mu0, no_ts, T, rec);
  // Print table
  o << std::setw(16) << "no. t.s." << std::setw(16) << "H1 seminorm"
    << std::setw(16) << "L2 norm bd\n";
  for (int j = 0; j < it_norms.size(); ++j) {
    auto [l2b, h1s] = it_norms[j];
    o << j << std::setw(16) << h1s << std::setw(16) << l2b << " ; ... \n";
  }
  return it_norms;
}
/* SAM_LISTING_END_6 */

}  // namespace IRKDegenerateEvl

#endif
