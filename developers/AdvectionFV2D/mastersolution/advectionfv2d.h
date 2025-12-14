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

#if SOLUTION
  int bound_nnz = dofh.Mesh()->NumEntities(0) * 4;
  B_Matrix.reserve(bound_nnz);

  // Iterate over all cells
  for (const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) {
    // Compute area of cell
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    double area = lf::geometry::Volume(*geo_p);

    // Corresponting DOF of ref. cell
    int row = dofh.GlobalDofIndices(*cell)[0];

    // Corresponting normal vectors of ref. cell
    Eigen::Matrix<double, 2, Eigen::Dynamic> cur_normal_vectors =
        (*normal_vectors)(*cell);

    // Get adjacent cells of ref. cell
    std::array<const lf::mesh::Entity*, 4> neighbour_cells =
        (*adjacentCells)(*cell);

    // Get edges of ref. cell
    auto cell_edges = cell->SubEntities(1);

    // Iterate over all edges
    int num_edges = cur_normal_vectors.cols();
    for (int edge_nr = 0; edge_nr < num_edges; ++edge_nr) {
      // Geo pointer of current edge
      const lf::geometry::Geometry* edge_geo_p =
          (cell_edges[edge_nr])->Geometry();

      // Length of current edge
      double edge_length = lf::geometry::Volume(*edge_geo_p);

      // Corners of current edge
      Eigen::Matrix<double, 2, Eigen::Dynamic> corner_edges =
          lf::geometry::Corners(*edge_geo_p);

      // Midpoint of current edge
      Eigen::Vector2d midpoint =
          (corner_edges.col(0) + corner_edges.col(1)) * 0.5;

      // Get neighbor cell
      const lf::mesh::Entity* next_cell = neighbour_cells[edge_nr];

      if (next_cell != nullptr) {
        // Corresponting DOF of next_cell
        int col = dofh.GlobalDofIndices(*next_cell)[0];

        // Compute Flux and store it directly in B
        double flux = (cur_normal_vectors.col(edge_nr)).dot(beta(midpoint));
        if (flux >= 0) {
          B_Matrix.coeffRef(row, row) -= flux * edge_length / area;
        } else {
          B_Matrix.coeffRef(row, col) -= flux * edge_length / area;
        }
      } else {
        // In case a cell has no neighbor at the current edge (boundary elem.),
        // the flux has to be considered if it is greater than 0

        // Compute Flux and store it directly in B
        double flux = cur_normal_vectors.col(edge_nr).dot(beta(midpoint));
        if (flux >= 0) {
          B_Matrix.coeffRef(row, row) -= flux * edge_length / area;
        }
      }
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

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

#if SOLUTION
  // Compute B
  Eigen::SparseMatrix B_matrix =
      initializeMOLODEMatrix(dofh, beta, adjacentCells, normal_vectors);

  // Enforce dirichlet boundary condition
  // Might interfere with the initial bump!
  // Flag Edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_blags{
      lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 1)};

  // Iterate over all cells
  for (const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) {
    int counter = 0;
    // Iterate over all edges
    for (auto edge : cell->SubEntities(1)) {
      // Cell contains edge at boundary
      if (bd_blags(*edge)) {
        // Corners of current edge
        Eigen::Matrix<double, 2, Eigen::Dynamic> corner_edges =
            lf::geometry::Corners(*edge->Geometry());

        // Midpoint of current edge
        Eigen::Vector2d midpoint =
            (corner_edges.col(0) + corner_edges.col(1)) * 0.5;

        // Corresponting normal vectors of ref. cell
        Eigen::Matrix<double, 2, Eigen::Dynamic> cur_normal_vectors =
            (*normal_vectors)(*cell);

        // If the flux is < 0, set entity in mu and the row of M to zero
        if ((cur_normal_vectors.col(counter)).dot(beta(midpoint)) < 0) {
          int idx = dofh.GlobalDofIndices(*cell)[0];
          B_matrix.row(idx) *= 0;
          mu[idx] = 0.0;
        }
      }
    }
    counter++;
  }

  // Timestepping
  Eigen::VectorXd k0;
  Eigen::VectorXd k1;
  const double overflow = 1000 * u0_h.lpNorm<Eigen::Infinity>();
  for (int step = 0; step < M; ++step) {
    double tau = T / M;
    k0 = B_matrix * mu;
    k1 = B_matrix * (mu + tau * k0);
    mu = mu + tau * 0.5 * (k0 + k1);

    if (mu.lpNorm<Eigen::Infinity>() > overflow) {
      throw std::overflow_error("Overflow occured!!\n");
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

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

#if SOLUTION
  // Inititialize a vector for the initial condition
  Eigen::VectorXd u0_h(num_dof);

  // Iterate over all cells
  for (const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Compute barycenter of the cell
    Eigen::Vector2d x = barycenter(corners);

    // Find the index of the DOF for the cell
    int idx = dofh.GlobalDofIndices(*cell)[0];

    u0_h[idx] = u0(x);
  }

  // Set up the number of steps according to the CFL condition
  int M = int((T / computeHmin(dofh.Mesh())) + 2);

  return solveAdvection2D(dofh, beta, u0_h, adjacentCells, normal_vectors, T,
                          M);
#else
  //====================
  // Your code goes here
  return Eigen::VectorXd::Zero(num_dof);
  //====================
#endif
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

#if SOLUTION
  // Iterate over all cells
  for (const lf::mesh::Entity* cell : dofh.Mesh()->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Compute the barycenter of the cell
    Eigen::Vector2d x = barycenter(corners);

    // Find the index of the DOF for the cell
    int idx = dofh.GlobalDofIndices(*cell)[0];

    // Compute exact solution
    Eigen::Vector2d phi_inv_x = phi_inv * x;
    ref_solution[idx] = u0(phi_inv_x);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

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
#if SOLUTION
  // Set upper and lower limit
  int M_upper = int((T / computeHmin(dofh.Mesh())) + 1);
  int M_lower = 1;

  // Compute cell normals
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      normal_vectors = AdvectionFV2D::computeCellNormals(dofh.Mesh());

  // Compute adjecent cells
  std::shared_ptr<
      lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(dofh.Mesh());

  // Initialize a vector for the result and
  // randomly initialize a vector for the initial condition
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd u0_h = Eigen::VectorXd::Random(num_dof);

  // Shift upper and lower bound until contition below isn't satisfied
  while ((M_upper - M_lower) > 2) {
    // Perform a simulation at M_middle
    // If it succeeds, set the upper bound to M_middle
    // Otherwise set the lower bound to M_middle
    // Repeat ...
    int M_middle = (M_upper + M_lower) / 2;
    try {
      solveAdvection2D(dofh, beta, u0_h, adjacentCells, normal_vectors, T,
                       M_middle);
      M_upper = M_middle;
    } catch (const std::exception& e) {
      M_lower = M_middle;
    }
  }

  int thres = M_upper;
  return thres;
#else
  //====================
  // Your code goes here
  // Replace the dummy return value below:
  return 0;
  //====================
#endif
}
/* SAM_LISTING_END_5 */

}  // namespace AdvectionFV2D

#endif  // define ADVECTIONFV2D_H
