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
    //====================
    // Your code goes here
    //====================
  }
  /* Public member functions */
  void compTimestep(double tau, Eigen::VectorXd& p, Eigen::VectorXd& q) const;
  double computeEnergies(const Eigen::VectorXd& p,
                         const Eigen::VectorXd& q) const;

 private:
  //====================
  // Your code goes here
  //====================
};  // class SympTimestepWaveEq
/* SAM_LISTING_END_3 */

/* Implementing member functions of class SympTimestepWaveEq */
/* SAM_LISTING_BEGIN_9 */
template <typename FUNCTION>
void SympTimestepWaveEq<FUNCTION>::compTimestep(double tau, Eigen::VectorXd& p,
                                                Eigen::VectorXd& q) const {
  //====================
  // Your code goes here
  //====================
}  // SympTimestepWaveEq<FUNCTION>::compTimestep
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_0 */
template <typename FUNCTION>
double SympTimestepWaveEq<FUNCTION>::computeEnergies(
    const Eigen::VectorXd& p, const Eigen::VectorXd& q) const {
  double energy;
  //====================
  // Your code goes here
  //====================
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
  //====================
  // Your code goes here
  //====================
  return solution_pair;
}
/* SAM_LISTING_END_7 */

void wavePropSimulation(unsigned int m);

double testStab();

}  // namespace SymplecticTimesteppingWaves

#endif
