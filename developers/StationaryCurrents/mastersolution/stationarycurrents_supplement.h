/** @file
 * @brief Extra definitions/function for StationaryCurrents
 * @author Ralf Hiptmair
 * @date August 2020
 * @copyright MIT License
 */

#ifndef DMXBCS_H_
#define DMXBCS_H_

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
#if SOLUTION
/** @brief Computation of gradients of barycentric coordinate functions,
 * exterior edge-length-weighted normals , and of the area of a triangle.
 *
 * @param vertices 2x3 matrix whose columns contain the vertex coordinates of
 * the triangle
 * @return tuple containing the gradients of barycentric coordinate functions,
 * exterior edge-length-weighted normals , and of the area of a triangle. The
 * gradients and the normals are stored in the columns of a 2x3 matrix.
 */
std::tuple<Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 2, 3>, double>
getTriangleGradLambdaNormals(Eigen::Matrix<double, 2, 3> vertices);

/** @brief Compute exterior edge-length-weighted normals for triangular or
 * quadrilateral cells
 *
 * @param corners 2xn-matrix whose columns contain the vertex coordinates of the
 * n vertices
 */
Eigen::MatrixXd exteriorCellNormals(const Eigen::MatrixXd& corners);

/** @brief Compute exterior edge-length-weighted normals for edges on the
 * boundary
 *
 * @param mesh_p pointer to finite element mesh containing only triangles with
 * straight edges.
 * @return Edge-indeed array containing edge-length-weighted normals for every
 * edge on the boundary , the zero vector for internal edges
 */
lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> exteriorEdgeWeightedNormals(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

// A debugging function
bool validateNormals(const lf::mesh::Mesh& mesh);
void printNodeTags(const lf::mesh::Mesh& mesh,
                   lf::mesh::utils::CodimMeshDataSet<int>& nodeids);

/** @brief Evaluation of boundary formula for contact flux
 *
 * @tparam SIGMAFUNCTION functor type for diffusion coefficient
 * @param fe_space pointer to lowest-order FE space object
 * @param sigma diffusion coefficient
 * @param edgeids edge-indexed array of id numbers marking contacts
 * @param contact_id
 *
 * Direct computation of flux through a contact based on a finite-element
 * solution of the mixed BVP.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename SIGMAFUNCTION>
double contactFlux(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma,
    const lf::mesh::utils::CodimMeshDataSet<int>& edgeids, int contact_id = 0) {
  // Obtain object managing dof indexing
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Check whether coefficient vector matches dof handler
  LF_ASSERT_MSG(sol_vec.size() == dofh.NumDofs(),
                "Size mismatch for coefficient vector");
  // Summation variable for returning result
  double s = 0.0;
#if SOLUTION
  // Counter for edges on selected contact (optional)
  unsigned int ed_cnt = 0;
#endif
  const lf::mesh::Mesh& mesh{*dofh.Mesh()};
  // We cannot loop over edges, because information from dofs not located on the
  // boundary is also required. Therefore we have to loop over all cells and
  // check whether they abut the relevant boundary part.
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Implemented for triangles only
    // Make sure the cell is of triangular shape
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "contactFlux: implemented for triangles only");
    // Obtain array of edge pointers (sub-entities of co-dimension 1)
    std::span<const lf::mesh::Entity* const> sub_ent_range{
        cell->SubEntities(1)};
    // Must be three edges
    LF_ASSERT_MSG(sub_ent_range.size() == 3, "Triangle must have three edges!");
    // Check whether a relevant contact edge belongs to the cell
    std::array<bool, 3> on_contact{false};
    unsigned int cnt = 0;
    // loop over the edges and check whether they belong to the boundary
    for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
      const lf::mesh::Entity& edge{*sub_ent_range[j]};
      if (edgeids(edge) == contact_id) {
        on_contact[j] = true;
        cnt++;
#if SOLUTION
        ed_cnt++;
#endif
      }
    }
    // The counter cnt counts the number of edges belonging to the selected
    // contact.
    if (cnt > 0) {
      // A contact edge belongs to the current cell
      // Compute the gradients, edge-weighted exterior normals and area
      const lf::geometry::Geometry& geo{*(cell->Geometry())};
      auto [grad_bary_coords, normals, area] =
          getTriangleGradLambdaNormals(lf::geometry::Corners(geo));
#if SOLUTION
      // Compute (constant) local gradient of the finite element solution
      // DofHandler must provide three degrees of freedom per cell for piecewise
      // linear Lagrangian finite elements on triangles
      LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                    "contactFlux: 3 dofs per triangle mandatory!");
      // Fetch an array of global indices of local shape functions
      const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
      // This is the gradient of the FE solution on the current cell
      const Eigen::Vector2d local_gradient{
          grad_bary_coords * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                              sol_vec[glob_dof_idx[1]],
                              sol_vec[glob_dof_idx[2]])
                                 .finished()};
      // Use midpoint local quadrature rule
      // Midpoints of edges in reference coordinates
      const Eigen::MatrixXd mp_ref{
          (Eigen::MatrixXd(2, 3) << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5).finished()};
      // Physical coordinates of midpoints of edges
      const Eigen::MatrixXd mp_phys{geo.Global(mp_ref)};
      // Sum contributions of local quadrature points
      for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
        if (on_contact[j]) {
          // (sigma*grad u_h)*normal
          s += normals.col(j).dot(sigma(mp_phys.col(j)) * local_gradient);
        }
      }
#else
      //====================
      // Your code goes here
      std::cerr << "SOLUTION MISSING" << std::endl;
      //====================
#endif
    }
  }  // end loop over cells
#if SOLUTION
  std::cout << "Summed flux for " << ed_cnt << " edges." << std::endl;
#endif
  return s;
}  // end contact flux
/* SAM_LISTING_END_2 */
#endif

#if SOLUTION
/** @brief Volume based formula for the evaluation of contact fluxes
 *
 * @param fe_space pointer to an object describing a lowest-order Lagrangian
 * finite element space
 * @param sol_vec basis expansion coefficient vector of FE solution
 * @param sigma diffusion coefficient: involved in the definition of the flux
 * @param gradpsi gradient of weighting function for volume formula.
 *
 * This function uses a quadrature which is exact for polynomials of degree 2
 */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFluxEXT(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  // Underlying FE mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Local-to-Global map for local/global shape function indices
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Obtain quadrature rule
  const lf::quad::QuadRule quadrule{
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2)};
  // Summation variable
  double s = 0.0;
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Check matching of reference element (unit triangle)
    LF_VERIFY_MSG(cell->RefEl() == quadrule.RefEl(),
                  "Mismatch of reference element for " << *cell);
    // Obtain geometry information for entity
    const lf::geometry::Geometry& geo{*cell->Geometry()};
    // Compute the gradients, edge-weighted exterior normals and area
    auto [grad_bary_coords, normals, area] =
        getTriangleGradLambdaNormals(lf::geometry::Corners(geo));
    // Compute (constant) local gradient of the finite element solution
    // DofHandler must provide three degrees of freedom per cell  piecewise
    // linear Lagrangian finite elements on triangles
    LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                  "contactFlux: 3 dofs per triangle mandatory!");
    const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
    const Eigen::Vector2d local_gradient{
        grad_bary_coords * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                            sol_vec[glob_dof_idx[1]], sol_vec[glob_dof_idx[2]])
                               .finished()};

    // Number of quadrature points
    const int P = quadrule.NumPoints();
    // Quadrature points
    const Eigen::MatrixXd zeta_ref{quadrule.Points()};
    // Map quadrature points to physical/world coordinates
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    // Quadrature weights
    const Eigen::VectorXd w_ref{quadrule.Weights()};
    // Gramian determinants
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};
    // Iterate over the quadrature points
    for (int l = 0; l < P; ++l) {
      const auto quadnode{zeta.col(l)};
      s += w_ref[l] * gram_dets[l] *
           (gradpsi(quadnode)).dot(sigma(quadnode) * local_gradient);
    }
  }
  return s;
}
#endif

#if SOLUTION
/** @see @ref stabFlux
 *
 * Alternative implementation based on @ref MeshFunction
 */
/* SAM_LISTING_BEGIN_5 */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFluxMF(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{fe_space->Mesh()};
  // Coefficient function and weight function
  const lf::mesh::utils::MeshFunctionGlobal mf_sigma(sigma);
  const lf::mesh::utils::MeshFunctionGlobal mf_gradpsi(gradpsi);
  // Build a MeshFunction representing the gradient of the finite element
  // solution
  const lf::fe::MeshFunctionGradFE mf_grad(fe_space, sol_vec);
  // Mesh function representing the integrand
  const auto mf_itg{lf::mesh::utils::transpose(mf_sigma * mf_grad) *
                    mf_gradpsi};
  const double s = lf::fe::IntegrateMeshFunction(
      *mesh_p, mf_itg, [](const lf::mesh::Entity& e) {
        return lf::quad::make_QuadRule(e.RefEl(), 2);
      })(0, 0);
  return s;
}  // end stabFluxMF
/* SAM_LISTING_END_5 */

/** @see @ref stabFlux
 *
 * Implementation based on midpoint quadrature rule
 */
/* SAM_LISTING_BEGIN_4 */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFluxMPR(
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
    // Compute the gradients, edge-weighted exterior normals and area
    // (The normals are not used here)
    // An alternative implementation could use ScalarReferenceElement
    auto [grad_bary_coords, normals, area] =
        getTriangleGradLambdaNormals(lf::geometry::Corners(geo));
    // DofHandler must provide three degrees of freedom per cell for piecewise
    // linear Lagrangian finite elements on triangles
    LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                  "contactFlux: 3 dofs per triangle mandatory!");
    // Fetch indices of global shape functions associated with triangle
    const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
    // Compute (constant) local gradient of the finite element solution
    const Eigen::Vector2d local_gradient{
        grad_bary_coords * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                            sol_vec[glob_dof_idx[1]], sol_vec[glob_dof_idx[2]])
                               .finished()};
    // Find physical/world coordinates of "midpoint" of triangle
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    // Add contribution of triangle
    const auto quadnode{zeta.col(0)};
    s += area * (gradpsi(quadnode)).dot(sigma(quadnode) * local_gradient);
  }
  return s;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_H */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFluxTRF(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  // Underlying FE mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Local-to-Global map for local/global shape function indices
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Reference coordinates of "midpoint" of a triangle (center of gravity)
  const Eigen::MatrixXd zeta_ref{
      (Eigen::Matrix<double, 2, 1>() << 1.0 / 3.0, 1.0 / 3.0).finished()};
  // Obtain gradients of reference shape functions at center of gravity
  const lf::fe::ScalarReferenceFiniteElement<double>& ref_lsf{
      *fe_space->ShapeFunctionLayout(lf::base::RefEl::kTria())};
  LF_ASSERT_MSG(
      ref_lsf.NumRefShapeFunctions() == 3,
      "There must be 3 reference shape functions associated with vertices");
  const Eigen::MatrixXd ref_lsf_grads{
      ref_lsf.GradientsReferenceShapeFunctions(zeta_ref)};
  LF_ASSERT_MSG(ref_lsf_grads.rows() == 3, "Three gradients required!");
  LF_ASSERT_MSG(ref_lsf_grads.cols() == 2, "Two components expected!");
  // Summation variable
  double s = 0.0;
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Not implemented for " << *cell);
    // Obtain geometry information for entity
    const lf::geometry::Geometry& geo{*cell->Geometry()};
    // Fetch the transformation matrix for gradients
    const Eigen::MatrixXd JinvT{geo.JacobianInverseGramian(zeta_ref)};
    LF_ASSERT_MSG(
        (JinvT.rows() == 2) && (JinvT.cols() == 2),
        "JinvT is " << JinvT.rows() << " x " << JinvT.cols() << "-matrix!");
    // Compute the gradients of the local shape functions by transformation
    const auto grad_lsf{JinvT * ref_lsf_grads.transpose()};
    // Compute the area of the triangle
    const double area = lf::geometry::Volume(geo);
    // DofHandler must provide three degrees of freedom per cell for piecewise
    // linear Lagrangian finite elements on triangles
    LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                  "contactFlux: 3 dofs per triangle mandatory!");
    // Fetch indices of global shape functions associated with triangle
    const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
    // Compute (constant) local gradient of the finite element solution
    const Eigen::Vector2d local_gradient{
        grad_lsf * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                    sol_vec[glob_dof_idx[1]], sol_vec[glob_dof_idx[2]])
                       .finished()};
    // Find physical/world coordinates of "midpoint" of triangle
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    // Add contribution of triangle
    const auto quadnode{zeta.col(0)};
    s += area * (gradpsi(quadnode)).dot(sigma(quadnode) * local_gradient);
  }
  return s;
}
/* SAM_LISTING_END_G */
#endif

}  // namespace dmxbc

#endif
