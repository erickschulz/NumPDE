/** @file
 * @brief Solving a second-order elliptic BVP with special mixed boundary
 * conditions
 * @author Ralf Hiptmair
 * @date July 2020
 * @copyright MIT License
 */

#ifndef DMXBC_H_
#define DMXBC_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>


namespace dmxbc {

/** @brief Computation of gradients of barycentric coordinate functions
 *
 * @param vertices 2x3 matrix whose columns contain the vertex coordinates of
 * the triangle
 * @return gradients of barycentric coordinate functions stored in the columns
 * of a 2x3 matrix.
 *
 */
Eigen::Matrix<double, 2, 3> GradsBaryCoords(
    Eigen::Matrix<double, 2, 3> vertices);

/** @brief Reading mesh from Gmsh file and flagging edges on contact boundary
 *         parts
 *
 * @param filename path to Gmsh mesh file (.msh-file)
 * @return pointer to read mesh and an array of integers index by the edges of
 * the mesh
 *
 * The .msh file must contain boundary edges with the two physical names
 * "Contact0" and "Contact1". In the returned array edges tagged with "Context0"
 * recive the id 0, edges belonging to the group "Contact1" get id 1. All other
 * edges are assigned -1.
 */
std::pair<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<int>>
readMeshWithTags(std::string filename);

/** @brief Spreading positive edge ids to endpoint nodes
 *
 * @param mesh_p pointer to underlying FE mesh
 * @param edgeids ids for the edges of the mesh
 * @return Array index by the nodes of the mesh containing their indices
 *
 * If an edge carries a positive id, then this id is set for both its endpoints.
 */
lf::mesh::utils::CodimMeshDataSet<int> tagNodes(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<int> edgeids);

/** @brief Solving mixed BVP based on lowest-order Lagrangian finite elements
 *
 * @tparam SIGMAFUNCTOR type for functor providing diffusion coefficient.
 * @param fe_space pointer to finite element space object
 * @param nodeflags array of consecutive ids 0,...,m for nodes of the mesh
 * @param voltvals array of values to which the finite element solution should
 * be set in the nodes with the corresponding id, which serves as an index into
 * the array.
 * @sigma diffusion coefficient
 *
 * Solution of a second-order elliptic scalar boundary value problem with fixed
 * solution values in particular nodes. Standard LehrFEM++ implementation.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename SIGMAFUNCTOR>
Eigen::VectorXd solveMixedBVP(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    lf::mesh::utils::CodimMeshDataSet<int>& nodeflags,
    std::vector<double>&& voltvals, SIGMAFUNCTOR&& sigma) {
  // Mesh functions for coefficients
  lf::mesh::utils::MeshFunctionGlobal mf_sigma{sigma};
  lf::mesh::utils::MeshFunctionConstant<double> mf_gamma{0.0};
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space = number of nodes of the mesh
  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_ASSERT_MSG(N_dofs == mesh.NumEntities(2),
                " N_dofs must agree with number of nodes");
  std::cout << "Solving BVP with " << N_dofs << " FE dofs" << std::endl;
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right-hand side vector
  Eigen::VectorXd phi(N_dofs);
  // Assembly of Galerkin matrix and right-hand side vector plus
  // special treatment of dofs at contacts
  //====================
  // Your code goes here
  //====================
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);  // LU decomposition
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  return sol_vec;
}  // end solveMixedBVP

/* SAM_LISTING_END_1 */

/** @see \ref contactFlux
 *
 * Alternative implementation making use of @ref MeshFunction
 */
/* SAM_LISTING_BEGIN_3 */
template <typename SIGMAFUNCTION>
double contactFluxMF(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma,
    const lf::mesh::utils::CodimMeshDataSet<int>& edgeids, int contact_id = 0) {
  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Variable for summing boundary flux
  double s = 0.0;
  // Counter for edges on selected contact
  unsigned int ed_cnt = 0;
//====================
// Your code goes here
//====================
  std::cout << "Summed flux for " << ed_cnt << " edges." << std::endl;
  return s;
}  // end contactFluxMF

/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_G */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFlux(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  // Underlying FE mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Local-to-Global map for local/global shape function indices
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Reference coordinates of "midpoint" of a triangle
  const Eigen::MatrixXd zeta_ref{
      (Eigen::Matrix<double, 2, 1>() << 1.0 / 3.0, 1.0 / 3.0).finished()};
  // Summation variable
  double s = 0.0;
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Not implemented for " << *cell);
    // Obtain geometry information for entity
    const lf::geometry::Geometry& geo{*cell->Geometry()};
    //====================
    // Your code goes here
    //====================
  }
  return s;
}
/* SAM_LISTING_END_G */

/** @brief Main driver function
 *
 * @param basename Filename of .msh-file without suffix
 */
std::tuple<double, double, double> computePotential(std::string basename);

}  // namespace dmxbc

#endif
