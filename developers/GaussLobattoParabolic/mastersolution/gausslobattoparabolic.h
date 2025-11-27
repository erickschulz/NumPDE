/**
 * @file gausslobattoparabolic.h
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef GAUSSLOBATTOPARABOLIC_H_
#define GAUSSLOBATTOPARABOLIC_H_

#include <lf/assemble/assemble.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <functional>

namespace GaussLobattoParabolic {

/**
 * @brief Compute the matrix \Blue{$\tilde{M}$} described in the exercise.
 *
 * @param fe_space finite element space, say of dimension \Blue{$N$}.
 * @return matrix of size \Blue{$N\times N$}
 */
lf::assemble::COOMatrix<double> initMbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space);

/**
 * @brief Compute the matrix \Blue{$\tilde{A}$} described in the exercise.
 *
 * @param fe_space finite element space, say of dimension \Blue{$N$}.
 * @return matrix of size \Blue{$N\times N$}
 */
lf::assemble::COOMatrix<double> initAbig(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space);

/* SAM_LISTING_BEGIN_3 */
class RHSProvider {
 public:
  // Disabled constructors
  RHSProvider() = delete;
  RHSProvider(const RHSProvider&) = delete;
  RHSProvider(RHSProvider&&) = delete;
  RHSProvider& operator=(const RHSProvider&) = delete;
  RHSProvider& operator=(const RHSProvider&&) = delete;
  // Main constructor; precomputations to be done here
  RHSProvider(const lf::assemble::DofHandler& dofh,
              std::function<double(double)> g);
  // Destructor
  virtual ~RHSProvider() = default;
  /* Class member functions */
  // Evaluation operator for right-hand-side vector
  Eigen::VectorXd operator()(double t) const;

 private:
  std::function<double(double)> g_;
#if SOLUTION
  Eigen::VectorXd zero_one_;
#else
  //====================
  // Other data members to be introduced by you
  //====================
#endif
};
/* SAM_LISTING_END_3 */

/**
 * @brief Compute the expansion coefficient vector \Blue{$\mu$} of the
 * solution at time T, corresponding to the basis of fe_space,
 * with vanishing inital data.
 *
 * @param fe_space finite element space, say of dimension \Blue{$N$}.
 * @param T final time
 * @param M number of time steps for the Gauss-Lobatto IIIC scheme
 * @param g continuously differentiable function satisfying \Blue{$g(0)=0$} that
 * describes the time-dependent boundary conditions.
 * @param u0 inital data satisfying \Blue{$g(0)=u_0$} on the boundary
 * @return vector of size \Blue{$N$} holding the coefficients of the solution
 */
/* SAM_LISTING_BEGIN_4 */
template <typename GFUNCTION>
Eigen::VectorXd evolveIBVPGaussLobatto(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    double T, unsigned int M, GFUNCTION&& g) {
  // timestep size
  const double tau = T / M;

  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  const int N = dofh.NumDofs();

  // Coefficient vector, initial value = 0
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(N);

  // Build left-hand side sparse block matrix
  lf::assemble::COOMatrix<double> lhs(2 * N, 2 * N);
#if SOLUTION
  // Arrange blocks by triplet manipulations
  // First: Two copies of \tilde{M} on the diagonal
  lf::assemble::COOMatrix<double> COO_M = initMbig(fe_space);
  for (const Eigen::Triplet<double>& triplet : COO_M.triplets()) {
    lhs.AddToEntry(triplet.row(), triplet.col(), triplet.value());
    lhs.AddToEntry(triplet.row() + N, triplet.col() + N, triplet.value());
  }
  // Scaled copies of \tilde{A} added to diagonal and set as off-diagonal blocks
  lf::assemble::COOMatrix<double> COO_A = initAbig(fe_space);
  for (const Eigen::Triplet<double>& triplet : COO_A.triplets()) {
    const int row = triplet.row();
    const int col = triplet.col();
    const double value = 0.5 * tau * triplet.value();
    lhs.AddToEntry(row, col, value);
    lhs.AddToEntry(row, col + N, -value);
    lhs.AddToEntry(row + N, col, value);
    lhs.AddToEntry(row + N, col + N, value);
  }
  // Use SparseLU for non-symmetric but square left-hand side matrix
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(lhs.makeSparse());

  // Get Sparse A to compute \Blue{$A\mu$} below
  Eigen::SparseMatrix<double> A = COO_A.makeSparse();
  // Helper class for computuing source term
  RHSProvider rhs_provider(dofh, std::forward<GFUNCTION>(g));

  // Main timesteping loop : initializations
  // Get time-dependent load vector for t=0
  Eigen::VectorXd phi = rhs_provider(0.0);

  // The actual timestepping loop
  for (unsigned int j = 0; j < M; ++j) {
    // Compute time-dependent right-hand side vector
    Eigen::VectorXd rhs(2 * N);
    Eigen::VectorXd A_mu = A * mu;
    rhs.head(N) = phi - A_mu;
    phi = rhs_provider(tau * (j + 1));
    rhs.tail(N) = phi - A_mu;

    // Compute increments
    Eigen::VectorXd k = solver.solve(rhs);

    // Perform a Gauss-Lobatto step to compute \Blue{$\vec{\mu}^{j+1}$}
    mu = mu + 0.5 * tau * (k.head(N) + k.tail(N));
  }
#else
  //====================
  // Your code goes here

  // Build the left-hand side matrix using initMbig(...) and initAbig(...)
  // replacing these dummy values:
  for (int i = 0; i < 2 * N; ++i) lhs.AddToEntry(i, i, 1.0);

  // Then ompute the LU decomposition outside of the timestepping loop:
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(lhs.makeSparse());
  // Helper class for computuing source term
  RHSProvider rhs_provider(dofh, std::move(g));

  // ...
  // Perform Timestepping here relying on the evaluation operators of
  // RHSProvider to obtain the time-dependent right-hand-side vector in
  // each timestep
  // ...
  //====================
#endif
  return mu;
}
/* SAM_LISTING_END_4 */

}  // namespace GaussLobattoParabolic

#endif  // #define GAUSSLOBATTOPARABOLIC_H_
