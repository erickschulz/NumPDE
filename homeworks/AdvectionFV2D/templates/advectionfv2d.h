#ifndef ADVECTIONFV2D_H
#define ADVECTIONFV2D_H
/**
 * @file advectionfv2d.h
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace AdvectionFV2D {

/**
 * @brief Compute the normal vectors for all edges of each cell.
 *
 * @param mesh_p Pointer to the mesh.
 * @return CodimMeshDataSet which stores for each cell all outer normal vectors.
 */
std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>>>
computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/**
 * @brief Find the neighbor cells to each cell in the mesh.
 *
 * @param mesh_p Pointer to the mesh.
 * @return CodimMeshDataSet which stores for each cell all neighbor cells.
 */
std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
getAdjacentCellPointers(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/**
 * @brief Get the coefficients of the barycentric coordinate functions for
 * a TRIA element.
 *
 * @param triangle Corners of the element.
 * @return Matrix providing the coefficients.
 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3>& triangle);

/**
 * @brief Compute the barycenter of a QUAD or TRIA
 *
 * @param corners Corners of the element.
 * @return Barycenter of the element.
 */
Eigen::Vector2d barycenter(const Eigen::MatrixXd corners);

/**
 * @brief Setup MOLODE Matrix
 *
 * @param dofh Reference to dof handler.
 * @param beta Functor to vector field beta
 * @param adjacentCells Pointer to CodimMeshDataSet containing information
 * about the neighbors of the cells
 * @param normal_vectors Pointer to CodimMeshdataSet containing information
 * about the normal vectors at the edges of the elements
 * @return MOLODE matrix.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename VECTORFIELD>
Eigen::SparseMatrix<double> initializeMOLODEMatrix(
    const lf::assemble::DofHandler& dofh, VECTORFIELD&& beta,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        std::array<const lf::mesh::Entity*, 4>>>
        adjacentCells,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors) {
  //////////////////////////////////////////////////////////////////////////////
  // TODO: Consistency Number of NNZ Entities
  //////////////////////////////////////////////////////////////////////////////
  // Set up matrix B
  int num_dof = dofh.NumDofs();
  Eigen::SparseMatrix<double> B_Matrix(num_dof, num_dof);

  //====================
  // Your code goes here
  //====================

  return B_Matrix;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute the minimum distance between the barycenters of two cells
 *
 * @param mesh_p mesh_p Pointer to the mesh.
 * @return Minimum distance.
 */
////////////////////////////////////////////////////////////////////////////////
// TODO: Inconsistency - Typo in computeMin arguments (constd -> const)
////////////////////////////////////////////////////////////////////////////////
double computeHmin(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

// Task 8-8.m
/**
 * @brief Function to initialize MOLODE Matrix, considering the boundary
 * condition and perform timestepping
 *
 * @param dofh Reference to dof handler.
 * @param beta Functor to vector field beta
 * @param u0_h Vector describing the initial bump
 * @param adjacentCells Pointer to CodimMeshDataSet containing information
 * about the neighbors of the cells
 * @param normal_vectors Pointer to CodimMeshdataSet containing information
 * about the normal vectors at the edges of the elements
 * @param T Final time
 * @param M Number of timesteps
 * @return Result after timestepping.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename VECTORFIELD>
Eigen::VectorXd solveAdvection2D(
    const lf::assemble::DofHandler& dofh, VECTORFIELD&& beta,
    const Eigen::VectorXd& u0_h,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        std::array<const lf::mesh::Entity*, 4>>>
        adjacentCells,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors,
    double T, unsigned int M) {
  // Set mu to inital bump
  Eigen::VectorXd mu = u0_h;

  //====================
  // Your code goes here
  //====================

  return mu;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Function computes the initial bump vector and the number of
 * required timesteps and also calls solveAdvection2D()
 *
 * @param dofh Reference to dof handler.
 * @param beta Functor to vector field beta
 * @param u0 Functor describing the initial bump
 * @param adjacentCells Pointer to CodimMeshDataSet containing information
 * about the neighbors of the cells
 * @param normal_vectors Pointer to CodimMeshdataSet containing information
 * about the normal vectors at the edges of the elements
 * @param T Final time
 * @return Result after timestepping.
 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR, typename VECTORFIELD>
Eigen::VectorXd simulateAdvection(
    const lf::assemble::DofHandler& dofh, VECTORFIELD&& beta, FUNCTOR&& u0,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        std::array<const lf::mesh::Entity*, 4>>>
        adjacentCells,
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors,
    double T) {
  // Get number of dofs
  int num_dof = dofh.NumDofs();

  //====================
  // Your code goes here
  return Eigen::VectorXd::Zero(num_dof);
  //====================
}
/* SAM_LISTING_END_3 */

/**
 * @brief Function returns the exact result of the specified problem
 *
 * @param dofh Reference to dof handler.
 * @param u0 Functor describing the initial bump
 * @param T Final time
 * @return Exact result.
 */
/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
Eigen::VectorXd refSolution(const lf::assemble::DofHandler& dofh, FUNCTOR&& u0,
                            double T) {
  // Setup inverted phi^-1
  Eigen::Matrix2d phi_inv;
  double t_mod = T / std::sqrt(2.0);
  phi_inv << std::cos(t_mod), std::sin(t_mod), -std::sin(t_mod),
      std::cos(t_mod);

  // Initialize vector for the result
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd ref_solution(num_dof);

  //====================
  // Your code goes here
  //====================

  return ref_solution;
}
/* SAM_LISTING_END_4 */

////////////////////////////////////////////////////////////////////////////////
// TODO: Returned int instead of double
////////////////////////////////////////////////////////////////////////////////
// Task 8-8.p
/**
 * @brief Function searching for the CFL threshold
 *
 * @param dofh Reference to dof handler.
 * @param beta Functor to vector field beta
 * @param T Final time
 * @return Minimum number of steps where no blowup occurs.
 */
/* SAM_LISTING_BEGIN_5 */
template <typename VECTORFIELD>
int findCFLthreshold(const lf::assemble::DofHandler& dofh, VECTORFIELD&& beta,
                     double T) {
  //====================
  // Your code goes here
  // Replace the dummy return value below:
  return 0;
  //====================
}
/* SAM_LISTING_END_5 */

}  // namespace AdvectionFV2D

#endif  // define ADVECTIONFV2D_H
