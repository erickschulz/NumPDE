#ifndef NPDECODES_RADAUTHREETIMESTEPPING_RADAUTHREETIMESTEPPING_H_
#define NPDECODES_RADAUTHREETIMESTEPPING_RADAUTHREETIMESTEPPING_H_

/** @file radauthreetimestepping.h
 * @brief NPDE RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

namespace RadauThreeTimestepping {

/**
 * @brief time depedent heat source
 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler& dofh,
                                    double time);

/**
 * @brief solve heat equation with rhsVectorHeat source as source
 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler& dofh,
                                   unsigned int m, double final_time);

/**
 * @brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin stiffness and mass matrices
 *
 * This function transforms every column and row
 * associated to a global index belonging to a degree of freedom lying on the
 * boundary to zero entries but the diagonal one which is set to 1.0
 *
 * @param selectvals The predicate identifying the boundary indices of the
 * rows and columns that are to be dropped
 */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsColumns(SELECTOR&& selectvals,
                           lf::assemble::COOMatrix<SCALAR>& A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  // Set the rows and columns of boundary DOFs to zero
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  // Set the diagonal entries of zeroed out rows and columns to 1
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/**
 * @brief This class implements a Lehrfem++ matrix provider defining a
 * LinFEMassMatrixProvider::Eval function returning the local MASS matrix for
 * linear first-order lagrange FE bases over triangular mesh (only!).
 * Integration over the triangular cells is performed using the trapezoidal
 * rule.
 */
class LinFEMassMatrixProvider {
 public:
  /**
   * @brief default constructor
   */
  explicit LinFEMassMatrixProvider() = default;

  /**
   * @brief Default implement: all cells are active
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief Main method for computing the element vector
   * @param cell refers to current cell for which the element vector is desired
   * @returns The mass matrix corresponding to the given cell
   *
   * The implementation uses an analytic formula defined over triangular cells
   **/
  Eigen::Matrix3d Eval(const lf::mesh::Entity& tria);
};

/**
 * @brief This class implements a Lehrfem++ matrix provider defining a
 * TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval function returning the local
 * contribution to the element vectors for linear first-order lagrange FE bases
 * over triangular mesh (only!). Integration over the triangular cells is
 * performed using the trapezoidal rule.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>  // lambda predicate
class TrapRuleLinFEElemVecProvider {
 public:
  /**
   * @brief Constructor storing the right hand side function
   */
  explicit TrapRuleLinFEElemVecProvider(FUNCTOR f) : f_(f) {}

  /**
   * @brief Default implement: all cells are active
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief Main method for computing the element vector
   * @param cell current cell for which the element vector is desired
   * @returns The element vector corresponding to the given cell
   *
   * The implementation uses simple vertex based quadrature and an approximation
   * of the volume of a cell just using the integration element at the
   * barycenter.*/
  Eigen::Vector3d Eval(const lf::mesh::Entity& tria);

 private:
  // f_ provides the evaluation of the source function at coordinates
  FUNCTOR f_;
};
/* SAM_LISTING_END_2 */

// Deduction guide for TrapRuleLinFEElemVecProvider
template <typename FUNCTOR>
TrapRuleLinFEElemVecProvider(FUNCTOR) -> TrapRuleLinFEElemVecProvider<FUNCTOR>;

// TrapRuleLinFEElemVecProvider
/* Implementing member function Eval of class TrapRuleLinFEElemVecProvider*/
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::Vector3d TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& tria) {
  Eigen::Vector3d ElemVec;
  //====================
  // Your code goes here
  //====================
  return ElemVec;
}
/* SAM_LISTING_END_3 */

/**
 * @brief class providing timestepping for heat equation
 */
class Radau3MOLTimestepper {
 public:
  // Disabled constructors
  Radau3MOLTimestepper() = delete;
  Radau3MOLTimestepper(const Radau3MOLTimestepper&) = delete;
  Radau3MOLTimestepper(Radau3MOLTimestepper&&) = delete;
  Radau3MOLTimestepper& operator=(const Radau3MOLTimestepper&) = delete;
  Radau3MOLTimestepper& operator=(const Radau3MOLTimestepper&&) = delete;

  // Main constructor; precomputations are done here
  Radau3MOLTimestepper(const lf::assemble::DofHandler& dofh);

  // Destructor
  virtual ~Radau3MOLTimestepper() = default;

  /* Class member functions */
  // Discrete evolution operator for Radau IIA 3rd order
  Eigen::VectorXd discreteEvolutionOperator(double time, double tau,
                                            const Eigen::VectorXd& mu) const;

 private:
  const lf::assemble::DofHandler& dofh_;  // dangerous
                                          //====================
                                          // Your code goes here
                                          //====================
};

}  // namespace RadauThreeTimestepping

#endif
