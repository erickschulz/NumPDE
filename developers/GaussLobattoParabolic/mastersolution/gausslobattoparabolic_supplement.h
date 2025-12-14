/**
 * @file gausslobattoparabolic_supplement.h
 * @brief special version of function for convergence studies
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef GAUSSLOBATTOPARABOLIC_S_H_
#define GAUSSLOBATTOPARABOLIC_S_H_

#include <lf/assemble/assemble.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <functional>

namespace GaussLobattoParabolic {
#if SOLUTION
template <
    typename RECORDER = std::function<void(double, const Eigen::VectorXd&)>>
Eigen::VectorXd evolveIBVPGaussLobatto_u0(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    double T, unsigned int M, std::function<double(double)> g,
    std::function<double(Eigen::Vector2d)> u0,
    RECORDER&& rec = [](double, const Eigen::VectorXd&) {}) {
  lf::mesh::utils::MeshFunctionGlobal u0_mf{u0};
  Eigen::VectorXd mu = lf::uscalfe::NodalProjection(*fe_space, u0_mf);

  double tau = T / M;
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(M + 1, 0.0, T);

  // Build left-hand side sparse block matrix
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  int N = dofh.NumDofs();
  lf::assemble::COOMatrix<double> lhs(2 * N, 2 * N);

  lf::assemble::COOMatrix<double> COO_M = initMbig(fe_space);
  for (const Eigen::Triplet<double>& triplet : COO_M.triplets()) {
    lhs.AddToEntry(triplet.row(), triplet.col(), triplet.value());
    lhs.AddToEntry(triplet.row() + N, triplet.col() + N, triplet.value());
  }

  lf::assemble::COOMatrix<double> COO_A = initAbig(fe_space);
  for (const Eigen::Triplet<double>& triplet : COO_A.triplets()) {
    int row = triplet.row(), col = triplet.col();
    double value = 0.5 * tau * triplet.value();
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

  // Get time-dependent load vector
  RHSProvider rhs_provider(dofh, std::move(g));
  Eigen::VectorXd phi = rhs_provider(t(0));
  rec(t(0), mu);

  for (unsigned int j = 0; j < M; ++j) {
    // Compute time-dependent right-hand side vector
    Eigen::VectorXd rhs(2 * N);
    Eigen::VectorXd A_mu = A * mu;
    rhs.head(N) = phi - A_mu;
    phi = rhs_provider(t(j + 1));
    rhs.tail(N) = phi - A_mu;

    // Compute increments
    Eigen::VectorXd k = solver.solve(rhs);

    // Perform a Gauss-Lobatto step to compute \Blue{$\vec{\mu}^{j+1}$}
    mu = mu + 0.5 * tau * (k.head(N) + k.tail(N));
    rec(t(j + 1), mu);
  }
  return mu;
}
#endif

}  // namespace GaussLobattoParabolic
