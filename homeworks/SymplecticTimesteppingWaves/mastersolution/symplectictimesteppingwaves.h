#ifndef SYMPLECTIC_HPP
#define SYMPLECTIC_HPP
/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

// homework includes
#include "symplectictimesteppingwaves_assemble.h"

namespace SymplecticTimesteppingWaves {

class progress_bar {
  static const auto overhead = sizeof " [100%]";
  std::ostream& os;
  const std::size_t bar_width;
  std::string message;
  const std::string full_bar;

 public:
  progress_bar(std::ostream& os, std::size_t line_width, std::string message_,
               const char symbol = '.')
      : os{os},
        bar_width{line_width - overhead},
        message{std::move(message_)},
        full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')} {
    if (message.size() + 1 >= bar_width ||
        message.find('\n') != std::string::npos) {
      os << message << '\n';
      message.clear();
    } else {
      message += ' ';
    }
    write(0.0);
  }

  progress_bar(const progress_bar&) = delete;
  progress_bar& operator=(const progress_bar&) = delete;

  ~progress_bar() {
    write(1.0);
    os << '\n';
  }

  void write(double fraction);
};

/** @brief This class precomputes the Galerkin matrices and the Eigen solver
 * based on the Cholesky decomposition (LDLT) that we use to perform symplectic
 * timestepping for the wave equaton (hyperbolic PDE) */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTION>
class SympTimestepWaveEq {
 public:
  /* Constructor */
  SympTimestepWaveEq(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
      FUNCTION c) {
    /* Creating the Galerkin Matrices for the wave equation*/
    // Assembling the element Galerkin matrix for the volume integrals
    auto one_coeff = [](Eigen::Vector2d x) -> double { return 1.0; };
    auto zero_coeff = [](Eigen::Vector2d x) -> double { return 0.0; };
    A_ = assembleGalerkinMatrix(fe_space_p, one_coeff, c, zero_coeff);
    // Assembling the Galerkin boundary mass matrix for the boundary integral
    // This requires a new call to assembleGalerkinMatrix using new coefficient
    // funtions so that only the mass integrals remains
    M_ = assembleGalerkinMatrix(fe_space_p, zero_coeff, one_coeff, zero_coeff);
    /* Precompute sparse Cholesky decomposition for the boundary mass matrix */
    solver_M_.compute(M_);
    LF_VERIFY_MSG(solver_M_.info() == Eigen::Success,
                  "Cholesky LDLT decomposition for sparse matrix M_ failed");
  }
  /* Public member functions */
  void compTimestep(double tau, Eigen::VectorXd& p, Eigen::VectorXd& q) const;
  double computeEnergies(const Eigen::VectorXd& p,
                         const Eigen::VectorXd& q) const;

 private:
  Eigen::SparseMatrix<double> A_;  // Galerkin matrix for volume integrals
  Eigen::SparseMatrix<double> M_;  // Galerkin Matrix for boundary integral
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_M_;
};  // class SympTimestepWaveEq
/* SAM_LISTING_END_3 */

/* Implementing member functions of class SympTimestepWaveEq */
/* SAM_LISTING_BEGIN_9 */
template <typename FUNCTION>
void SympTimestepWaveEq<FUNCTION>::compTimestep(double tau, Eigen::VectorXd& p,
                                                Eigen::VectorXd& q) const {
  // Coefficients of the method
  Eigen::VectorXd a(3);
  a << 2.0 / 3.0, -2.0 / 3.0, 1.0;
  Eigen::VectorXd b(3);
  b << 7.0 / 24.0, 3.0 / 4.0, -1.0 / 24.0;
  // Solving for f(q,t) using precomputed Galerkin matrices and the
  // stored factorization of the mass matrix
  Eigen::VectorXd fq_in = -solver_M_.solve(A_ * q);
  // one step method
  for (int i = 0; i < 3; ++i) {
    p += tau * b(i) * fq_in;
    q += tau * a(i) * p;
    fq_in = -solver_M_.solve(A_ * q);
  }
}  // SympTimestepWaveEq<FUNCTION>::compTimestep
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_0 */
template <typename FUNCTION>
double SympTimestepWaveEq<FUNCTION>::computeEnergies(
    const Eigen::VectorXd& p, const Eigen::VectorXd& q) const {
  double energy;
  energy = 0.5 * (p.dot(M_ * p) + q.dot(A_ * q));
  return energy;
}  // SympTimestepWaveEq<FUNCTION>::computeEnergies
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_7 */
template <typename FUNCTION>
std::pair<Eigen::VectorXd, Eigen::VectorXd> solvewave(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fes_p,
    FUNCTION c, const Eigen::VectorXd& u0_vec, const Eigen::VectorXd& v0_vec,
    double T, unsigned int m) {
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solution_pair;
  double tau = T / m;  // time step
  std::cout << "Solving with uniform step size tau = " << tau << std::endl;
  progress_bar progress{std::clog, 70u, "Timestepping"};
  double progress_pourcentage;
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());  // dim. of FE space size
  LF_VERIFY_MSG(u0_vec.size() == N_dofs, "Wrong size of initial conditions u0");
  LF_VERIFY_MSG(v0_vec.size() == N_dofs, "Wrong size of initial conditions");

  // Precomputing the required data for symplectic time stepping
  SympTimestepWaveEq<decltype(c)> timestepper(fes_p, c);

  /* Starting the evolution using the initial conditions */
  Eigen::VectorXd q = u0_vec;
  Eigen::VectorXd p = v0_vec;
  // Initialization to store energy at different times
  Eigen::VectorXd energies(m + 1);

  /* Iterating symplectic stepping */
  for (int i = 0; i < m; i++) {
    energies[i] = timestepper.computeEnergies(p, q);
    timestepper.compTimestep(tau, p, q);
    // \textbf{Throw an exception} in case of severe increase of  the total
    // energy, which is a conserved quantity for the exact evolution.
    if (energies[i] > 10.0 * energies[0]) {
      throw "ENERGY BLOW UP";
    }
    progress_pourcentage = ((double)i) / m * 100.0;
    progress.write(progress_pourcentage / 100.0);
  }

  energies[m] = timestepper.computeEnergies(p, q);
  solution_pair = std::make_pair(q, energies);
  return solution_pair;
}
/* SAM_LISTING_END_7 */

void wavePropSimulation(unsigned int m);

double testStab();

}  // namespace SymplecticTimesteppingWaves

#endif
