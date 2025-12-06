/**
 * @file sufem.h
 * @brief NPDE homework SUFEM code
 * @author R. Hiptmair
 * @date July 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef SUFEM_H_
#define SUFEM_H_

#include <cmath>
#include <iostream>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/base/lf_assert.h>
#include <lf/base/ref_el.h>
#include <lf/fe/fe.h>
#include <lf/fe/loc_comp_ellbvp.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_traits.h>
#include <lf/mesh/utils/special_entity_sets.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/fe_space_lagrange_o1.h>
#include <lf/uscalfe/fe_space_lagrange_o2.h>
#include <lf/uscalfe/uscalfe.h>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace SUFEM {
/** @brief Entity matrix provider for advection bilinear form and Lagrangian
   finite elements */

template <typename VELOCITY>
class AdvectionElementMatrixProvider {
  static_assert(lf::mesh::utils::MeshFunction<VELOCITY>);

 public:
  using ElemMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  AdvectionElementMatrixProvider(const AdvectionElementMatrixProvider&) =
      delete;
  AdvectionElementMatrixProvider(AdvectionElementMatrixProvider&&) noexcept =
      default;
  AdvectionElementMatrixProvider& operator=(
      const AdvectionElementMatrixProvider&) = delete;
  AdvectionElementMatrixProvider& operator=(AdvectionElementMatrixProvider&&) =
      delete;

  /**
   * @brief Constructor: cell-independent precomputations and auto quadrature
   * rule
   *
   * @param fe_space collection of specifications for scalar-valued parametric
   * reference elements
   * @param velo mesh function for the advection velocity vector field
   *
   * This constructor uses local quadature rules with double the degree of
   * exactness as the polynomial degree of the finite element space.
   */
  AdvectionElementMatrixProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      VELOCITY velo);
  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   *
   * @throw lf::base::LfException in case the finite element specification is
   * missing for the type of the cell or if there is no quadrature rule
   * specified for the given cell type.
   */
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity& cell);

  /** Virtual destructor */
  virtual ~AdvectionElementMatrixProvider() = default;

 private:
  /** @name functors providing coefficient functions
   * @{ */
  /** Advection velocity coefficient */
  VELOCITY velo_;
  /** @} */
  // fe_precomp_[i] contains precomputed reference finite element for ref_el i.
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>, 5>
      fe_precomp_;
};

// Constructor (automatic choice of quadrature rules)
template <typename VELOCITY>
AdvectionElementMatrixProvider<VELOCITY>::AdvectionElementMatrixProvider(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    VELOCITY velo)
    : velo_(std::move(velo)), fe_precomp_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    if (fe != nullptr) {
      fe_precomp_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
              fe, lf::quad::make_QuadRule(ref_el, 2 * fe->Degree()));
    }
  }
}

// Main method for the computation of the element matrix
/* SAM_LISTING_BEGIN_5 */
template <typename VELOCITY>
typename SUFEM::AdvectionElementMatrixProvider<VELOCITY>::ElemMat
AdvectionElementMatrixProvider<VELOCITY>::Eval(const lf::mesh::Entity& cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe =
      fe_precomp_[ref_el.Id()];
  if (!pfe.isInitialized()) {
    // Accident: cell is of a type not covered by finite element
    // specifications or there is no quadrature rule available for this
    // reference element type
    std::stringstream temp;
    temp << "No local shape function information or no quadrature rule for "
            "reference element type "
         << ref_el;
    throw lf::base::LfException(temp.str());
  }

  // Query the shape of the cell
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");

  // Physical dimension of the cell (must be 2)
  const lf::base::dim_t world_dim = geo_ptr->DimGlobal();
  LF_ASSERT_MSG_CONSTEXPR(world_dim == 2, "Only available for flat domains");
  // Gram determinant at quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      determinants.size() == pfe.Qr().NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << pfe.Qr().NumPoints());
  // Fetch the transformation matrices for the gradients
  const Eigen::MatrixXd JinvT(
      geo_ptr->JacobianInverseGramian(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      JinvT.cols() == 2 * pfe.Qr().NumPoints(),
      "Mismatch " << JinvT.cols() << " <-> " << 2 * pfe.Qr().NumPoints());
  LF_ASSERT_MSG(JinvT.rows() == world_dim,
                "Mismatch " << JinvT.rows() << " <-> " << world_dim);

  // Get the velocity vectors at quadrature points in the cell
  auto veloval = velo_(cell, pfe.Qr().Points());

  // Element matrix
  ElemMat mat(pfe.NumRefShapeFunctions(), pfe.NumRefShapeFunctions());
  mat.setZero();

  // Loop over quadrature points
  for (lf::base::size_type k = 0; k < pfe.Qr().NumPoints(); ++k) {
    const double w = pfe.Qr().Weights()[k] * determinants[k];
    // Matrix $\VG$ whose columns contain transformed gradients of all local
    // shape functions
    const auto trf_grad(
        JinvT.block(0, 2 * static_cast<Eigen::Index>(k), world_dim, 2) *
        pfe.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.rows(), 2)
            .transpose());
    // Row vector of transformed gradients multiplied with velocity vectors
    /*********************************************************
    Your implementation to fill in 'mat' goes here
    *********************************************************/
  }
  return mat;
}
/* SAM_LISTING_END_5 */

/** @brief Auxliary class providing a MeshFunction realizing the diffusion
   tensor in SU bilinear form */
template <typename VELOCITY>
class MeshFunctionDiffTensor {
  static_assert(lf::mesh::utils::MeshFunction<VELOCITY>);

 public:
  MeshFunctionDiffTensor(const MeshFunctionDiffTensor&) = default;
  MeshFunctionDiffTensor(MeshFunctionDiffTensor&&) noexcept = default;
  MeshFunctionDiffTensor& operator=(const MeshFunctionDiffTensor&) = delete;
  MeshFunctionDiffTensor& operator=(MeshFunctionDiffTensor&&) = delete;
  explicit MeshFunctionDiffTensor(VELOCITY velo) : velo_(std::move(velo)) {}
  virtual ~MeshFunctionDiffTensor() = default;
  // Local evaluation operator
  [[nodiscard]] std::vector<Eigen::Matrix2d> operator()(
      const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const;

 private:
  VELOCITY velo_;
};

/* SAM_LISTING_BEGIN_7 */
template <typename VELOCITY>
std::vector<Eigen::Matrix2d> MeshFunctionDiffTensor<VELOCITY>::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const {
  std::vector<Eigen::Matrix2d> ret;  // Element matrix
  const lf::geometry::Geometry& geo{*e.Geometry()};
  const double area = lf::geometry::Volume(geo);
  const Eigen::MatrixXd& corners_refc(e.RefEl().NodeCoords());
  const auto velo_cvals(velo_(e, corners_refc));  // $\Vv$ at cell corners
  const auto velovals(velo_(e, local));           // $\Vv$ in given points

  // Formula \prbeqref{eq:delta}
  /*********************************************************
  Adjust the following loop to fill in 'ret' correctly
  *********************************************************/
  for (auto& velovec : velovals) {
    ret.push_back(Eigen::Matrix2d());
  }
  return ret;
}
/* SAM_LISTING_END_7 */

/** @brief Assembly of full Galerkin matrix for streamline upwind bilinear form
 */
/* SAM_LISTING_BEGIN_6 */
template <typename VELOCITY>
lf::assemble::COOMatrix<double> buildSUGalerkinMatrix(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    VELOCITY velo) {
  static_assert(lf::mesh::utils::MeshFunction<VELOCITY>);
  // Create object of helper class providing diffusion tensor for SU bilinear
  // form
  MeshFunctionDiffTensor mf_diff(velo);
  // Local computations for advection bilinear form
  AdvectionElementMatrixProvider adv_elmat(fe_space, velo);
  // Local computations for diffusive part
  lf::fe::DiffusionElementMatrixProvider diff_elmat(fe_space, mf_diff);
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs = dofh.NumDofs();
  // Galerkin matrix in triplet format
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  // Invoke cell-oriented assembly in turns for both parts of the bilinear form
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, adv_elmat, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, diff_elmat, A_COO);
  return A_COO;
}
/* SAM_LISTING_END_6 */

/** Mark mesh nodes located on the (closed) inflow boundary */
template <typename VELOCITY>
lf::mesh::utils::CodimMeshDataSet<bool> flagNodesOnInflowBoundary(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p, VELOCITY velo) {
  static_assert(lf::mesh::utils::MeshFunction<VELOCITY>);
  // Array for flags
  lf::mesh::utils::CodimMeshDataSet<bool> nd_inflow_flags(mesh_p, 2, false);
  // Reference coordinates of center of gravity of a triangle
  const Eigen::MatrixXd c_hat = Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0);
  // Reference coordinates of midpoints of edges
  const Eigen::MatrixXd mp_hat =
      (Eigen::Matrix<double, 2, 3>() << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5)
          .finished();
  // Find edges (codim = 1) on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> ed_bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1));
  // Run through all cells of the mesh and determine
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    // Fetch geometry object for current cell
    const lf::geometry::Geometry& K_geo{*(cell->Geometry())};
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Only implemented for triangles");
    LF_ASSERT_MSG(K_geo.DimGlobal() == 2, "Mesh must be planar");
    // Obtain physical coordinates of barycenter of triangle
    const Eigen::Vector2d center{K_geo.Global(c_hat).col(0)};
    // Get velocity values in the midpoints of the edges
    auto velo_mp_vals = velo(*cell, mp_hat);
    // Retrieve pointers to all edges of the triangle
    std::span<const lf::mesh::Entity* const> edges{cell->SubEntities(1)};
    LF_ASSERT_MSG(edges.size() == 3, "Triangle must have three edges!");
    for (int k = 0; k < 3; ++k) {
      if (ed_bd_flags(*edges[k])) {
        const lf::geometry::Geometry& ed_geo{*(edges[k]->Geometry())};
        const Eigen::MatrixXd ed_pts{lf::geometry::Corners(ed_geo)};
        // Direction vector of the edge
        const Eigen::Vector2d dir = ed_pts.col(1) - ed_pts.col(0);
        // Rotate counterclockwise by 90 degrees
        const Eigen::Vector2d ed_normal = Eigen::Vector2d(dir(1), -dir(0));
        // For adjusting direction of normal so that it points into the exterior
        // of the domain
        const int ori = (ed_normal.dot(center - ed_pts.col(0)) > 0) ? -1 : 1;
        // Check angle of exterior normal and velocity vector
        const int v_rel_ori =
            ((velo_mp_vals[k].dot(ed_normal) > 0) ? 1 : -1) * ori;
        if (v_rel_ori < 0) {
          // Inflow: obtain endpoints of the edge and mark them
          std::span<const lf::mesh::Entity* const> endpoints{
              edges[k]->SubEntities(1)};
          LF_ASSERT_MSG(endpoints.size() == 2, "Edge must have two endpoints!");
          nd_inflow_flags(*endpoints[0]) = true;
          nd_inflow_flags(*endpoints[1]) = true;
        }
      }
    }
  }
  return nd_inflow_flags;
}

/** Solve Dirichlet problem for pure advection with no source */
/* SAM_LISTING_BEGIN_8 */
template <typename VELOCITY, typename DIRICHLET_DATA>
Eigen::VectorXd solveAdvectionDirichlet(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    VELOCITY velo, DIRICHLET_DATA mf_g) {
  static_assert(lf::mesh::utils::MeshFunction<VELOCITY>);
  static_assert(lf::mesh::utils::MeshFunction<DIRICHLET_DATA>);
  // Obtain full Galerkin matrix in COO format
  lf::assemble::COOMatrix<double> A_COO =
      SUFEM::buildSUGalerkinMatrix(fe_space, velo);
  // Zero right-hand side vector
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs = dofh.NumDofs();
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
  // ** Set boundary conditions **
  // Interpolate Dirichlet data
  Eigen::VectorXd g_coeffs = lf::fe::NodalProjection(*fe_space, mf_g);
  // Determine nodes on the inflow boundary part
  auto inflow_nodes{SUFEM::flagNodesOnInflowBoundary(fe_space->Mesh(), velo)};
  /*********************************************************
  Your implementation to adjust A_COO and phi to enforce Dirichlet BCs goes here
  *********************************************************/
  // Solve linear system using sparse elimination
  Eigen::SparseMatrix<double> A(A_COO.makeSparse());
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not factorize A");
  }
  return solver.solve(phi);
}
/* SAM_LISTING_END_8 */

void testSUFEMConvergence(unsigned int reflevels = 6,
                          const char* filename = nullptr);

}  // namespace SUFEM

#endif
