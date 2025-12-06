/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <unsupported/Eigen/KroneckerProduct>

namespace SDIRKMethodOfLines {

/** @brief class providing timestepping for convective cooling problem */
/* SAM_LISTING_BEGIN_1 */
class SDIRK2Timestepper {
 public:
  // Disabled constructor
  SDIRK2Timestepper() = delete;
  SDIRK2Timestepper(const SDIRK2Timestepper&) = delete;
  SDIRK2Timestepper(SDIRK2Timestepper&&) = delete;
  SDIRK2Timestepper& operator=(const SDIRK2Timestepper&) = delete;
  SDIRK2Timestepper& operator=(const SDIRK2Timestepper&&) = delete;
  // Main constructor; precomputations are done here
  explicit SDIRK2Timestepper(const lf::assemble::DofHandler& dofh, double tau,
                             double cool_coeff);
  // Destructor
  virtual ~SDIRK2Timestepper() = default;

  /* Class member functions */
  // Discrete evolution operator for SDIRK-2
  Eigen::VectorXd discreteEvolutionOperator(const Eigen::VectorXd& mu) const;

 private:
  double tau_;  // step size (in time)
  //====================
  // Your code goes here
  //====================
};
/* SAM_LISTING_END_1 */

/* Declaration of the functions of the library sdirkmethodoflines.h */
double thermalEnergy(const lf::assemble::DofHandler&, const Eigen::VectorXd&);

std::pair<Eigen::VectorXd, Eigen::VectorXd> solveTemperatureEvolution(
    const lf::assemble::DofHandler&, unsigned int, double, Eigen::VectorXd);

std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler& dofh,
                         double cool_coeff);

}  // namespace SDIRKMethodOfLines
