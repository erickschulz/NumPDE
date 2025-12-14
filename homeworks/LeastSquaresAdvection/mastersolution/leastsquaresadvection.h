/**
 * @file leastsquaresadvection.h
 * @brief NPDE homework LeastSquaresAdvection code
 * @author R. Hiptmair
 * @date July 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef LeastSquaresAdvection_H_
#define LeastSquaresAdvection_H_

// Include almost all parts of LehrFEM++; soem my not be needed
#include <lf/assemble/assemble.h>
#include <lf/assemble/assembly_types.h>
#include <lf/assemble/fix_dof.h>
#include <lf/base/lf_assert.h>
#include <lf/fe/fe.h>
#include <lf/fe/fe_tools.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/mesh/utils/mesh_function_unary.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace LeastSquaresAdvection {

/** @brief ENTITY MATRIX PROVIDER class for advection least squares bilinear
 * form
 *
 *
 */
template <lf::mesh::utils::MeshFunction REACTION_COEFF>
class LSQAdvectionMatrixProvider {
 public:
  using Scalar = double;
  using ElemMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  // Standard constructors
  LSQAdvectionMatrixProvider(const LSQAdvectionMatrixProvider&) = delete;
  LSQAdvectionMatrixProvider(LSQAdvectionMatrixProvider&&) noexcept = default;
  LSQAdvectionMatrixProvider& operator=(const LSQAdvectionMatrixProvider&) =
      delete;
  LSQAdvectionMatrixProvider& operator=(LSQAdvectionMatrixProvider&&) = delete;
  virtual ~LSQAdvectionMatrixProvider() = default;
  // Constructor, initializes data members and cell-indepedent quantities
  LSQAdvectionMatrixProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      Eigen::Vector2d velocity, REACTION_COEFF kappa);
  // Standard interface for ENTITY MATRIX PROVIDERS
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  ElemMat Eval(const lf::mesh::Entity& cell);

 private:
  REACTION_COEFF kappa_;      // Reaction coefficient $\cob{\kappa=\kappa(\Bx)}$
  Eigen::Vector2d velocity_;  // Velocity vector $\cob{\Vv}$
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>, 5>
      fe_precomp_;
};

template <lf::mesh::utils::MeshFunction REACTION_COEFF>
LSQAdvectionMatrixProvider<REACTION_COEFF>::LSQAdvectionMatrixProvider(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::Vector2d velocity, REACTION_COEFF kappa)
    : velocity_(velocity), kappa_(std::move(kappa)), fe_precomp_() {
  // Copied from loc_comp_ellbvp.h, constructor of
  // ReactionDiffusionElementMatrixProvider
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    // Check whether shape functions for that entity type are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe != nullptr) {
      // Precompute cell-independent quantities based on quadrature rules
      // with twice the degree of exactness compared to the degree of the
      // finite element space.
      fe_precomp_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
              fe, lf::quad::make_QuadRule(ref_el, 2 * fe->Degree()));
    }
  }
}

/* SAM_LISTING_BEGIN_1 */
template <lf::mesh::utils::MeshFunction REACTION_COEFF>
typename LSQAdvectionMatrixProvider<REACTION_COEFF>::ElemMat
LSQAdvectionMatrixProvider<REACTION_COEFF>::Eval(const lf::mesh::Entity& cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  const lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe =
      fe_precomp_[ref_el.Id()];
  LF_ASSERT_MSG(pfe.isInitialized(), "Precomputed data missing!");
  // Query the shape of the cell
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  // Physical dimension of the cell
  const lf::base::dim_t world_dim = geo_ptr->DimGlobal();
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
  // Element matrix
  ElemMat mat(pfe.NumRefShapeFunctions(), pfe.NumRefShapeFunctions());
  mat.setZero();
  // Compute values of $\cob{\kappa}$ in the quadrature nodes
  // $\cob{\Phibf_K(\wh{\zetabf}_{\ell})}$.
  auto kappaval = kappa_(cell, pfe.Qr().Points());
  // Loop over quadrature points
  for (lf::base::size_type k = 0; k < pfe.Qr().NumPoints(); ++k) {
    // The weighting factors $\cob{\nu_{\ell}}$
    const double w = pfe.Qr().Weights()[k] * determinants[k];
    // Transformed gradients in $\cob{\VG_{\ell}}$, see \prbeqref{eq:G}.
    const auto trf_grad(
        JinvT.block(0, 2 * static_cast<Eigen::Index>(k), world_dim, 2) *
        pfe.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.rows(), 2)
            .transpose());
    // Shape functions multiplied with kappa, vectors $\cob{\Bk_{\ell}}$,
    // \prbeqref{eq:kvec}
    const auto kvec = kappaval[k] * pfe.PrecompReferenceShapeFunctions().col(k);
    // Transformed gradients multiplied with velocity,
    // $\cob{\Vv^{\top}\VG_{\ell}}$
    const auto vTG(velocity_.transpose() * trf_grad);
    // Realization of \prbeqref{eq:Acomp}.
    const auto zvec = vTG + kvec.transpose();
    mat += w * zvec.transpose() * zvec;
  }
  return mat;
}
/* SAM_LISTING_END_1 */

/** @brief Mark mesh entities on inflow boundary
 *
 * @param mesh_p pointer to 2D hybrid finite element mesh
 * @param velocity constant velocity vector
 * @return Mesh data set with true flag for every mesh entitiy on/in the closure
 *         of the inflow boundary part
 */
lf::mesh::utils::AllCodimMeshDataSet<bool> flagEntitiesOnInflow(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    Eigen::Vector2d velocity);

/** @brief Compute least squares finite element Galerkin solution
 *
 * @tparam REACTION_COEFF mesh function type provoding reaction coefficient
 * @tparam GFUNCTION mesh function supplying Dirichlet data
 * @param fe_space pointer to underlying Lagrangian FE space
 * @param velocity constant velocity vector
 * @param kappa mesh function object for reaction coefficient
 * @param g mesh function object for boundary data
 */
/* SAM_LISTING_BEGIN_3 */
template <lf::mesh::utils::MeshFunction REACTION_COEFF,
          lf::mesh::utils::MeshFunction GFUNCTION>
Eigen::VectorXd solveAdvectionDirichletBVP(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::Vector2d velocity, const REACTION_COEFF& kappa, const GFUNCTION& g) {
  // I. Assemble the full Galerkin matrix for the least squares variational
  // formulation Fetch DofHandler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const size_t N = dofh.NumDofs();
  // Sparse matrix in triplet format
  lf::assemble::COOMatrix<double> A_coo(N, N);
  // Set up ENTITY\_MATRIX\_PROVIDER
  LeastSquaresAdvection::LSQAdvectionMatrixProvider lsq_adv_emp(
      fe_space, velocity, kappa);
  // Assembly of Galerkin matrix
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, lsq_adv_emp, A_coo);

  // II. Enforce Dirichlet boundary conditions on inflow boundary
  // as in \lref{sec:essbdc}
  const lf::mesh::utils::AllCodimMeshDataSet<bool> inflow_flags{
      flagEntitiesOnInflow(fe_space->Mesh(), velocity)};
  // Right-hand side vector
  Eigen::VectorXd phi(N);
  phi.setZero();
  // Interpolate extended boundary data to obtain $\cob{\wt{g}_h}$
  const Eigen::VectorXd gt_vec = lf::fe::NodalProjection(*fe_space, g);
  // Eliminate degrees of freedom on the boundary using
  // \lfpp's built-in function, see \lref{par:lfffsc}
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&inflow_flags, &dofh,
       &gt_vec](lf::assemble::glb_idx_t dof_idx) -> std::pair<bool, double> {
        const bool on_inflow = inflow_flags(dofh.Entity(dof_idx));
        if (on_inflow) {
          return {true, gt_vec[dof_idx]};
        }
        return {false, 0.0};
      },
      A_coo, phi);
  // III. Solving linear system of equations using
  // \eigen's sparse direct elimination solver
  Eigen::SparseMatrix<double> A_crs{A_coo.makeSparse()};
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd result = Eigen::VectorXd::Constant(N, 3.0);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("LU decomposition failed");
  }
  result = solver.solve(phi);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Solving LSE failed");
  }
  return result;
}
/* SAM_LISTING_END_3 */

/** @brief Tracking of L2 norm of Galerkin discretization error on sequence of
 * meshes obtained by regular refinement.
 *
 * @tparam LAGRFESPACE type derived from
 *          lf::uscalfe::UniformScalarFESpace<double>
 * @tparam functor type std::function<double(Eigen::Vector2d)>
 * @param g functor object providing Dirichlet data
 * @param kappa_val constant value of reaction coefficient
 * @param refsteps number of refinement steps to generate mesh hierarchy
 *
 * Function produces tabular output
 */

template <typename LAGRFESPACE, typename GFUNCTOR>
void testCVGLSQAdvectionReaction(GFUNCTOR&& g, double kappa_val,
                                 unsigned int refsteps = 4) {
  // Velocity vector (must not be changed!)
  const Eigen::Vector2d v(2.0, 1.0);
  //  ********** Part I: Manufactured solution  **********
  // Exact solution
  auto u_ex = [&](Eigen::Vector2d x) -> double {
    LF_ASSERT_MSG(
        (x[0] >= 0.0) and (x[0] <= 1.0) and (x[1] >= 0.0) and (x[1] <= 1.0),
        "x must lie inside the unit sqaure!");
    double tau;
    if (v[1] * x[0] > v[0] * x[1]) {
      tau = x[1] / v[1];
    } else {
      tau = x[0] / v[0];
    }
    Eigen::Vector2d x0 = x - tau * v;
    return g(x0) * std::exp(-kappa_val * tau);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_uex{u_ex};
  // Mesh function for reaction coefficient = 1
  lf::mesh::utils::MeshFunctionGlobal mf_kappa{
      [&kappa_val](Eigen::Vector2d /*x*/) -> double { return kappa_val; }};
  // Mesh function for boundary data
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};

  // ********** Part II: Loop over sequence of meshes **********
  // Generate a small unstructured triangular mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  const std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_ptr,
                                                              refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);

  // Number of levels
  const int L = multi_mesh.NumLevels();

  // Table of various error norms
  std::vector<std::pair<size_t, double>> errs;
  for (int level = 0; level < L; ++level) {
    const std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Piecewise quadrative Lagrangian finite element space
    auto fes_p = std::make_shared<LAGRFESPACE>(lev_mesh_p);
    // Solve on a single mesh
    Eigen::VectorXd mu_vec =
        solveAdvectionDirichletBVP(fes_p, v, mf_kappa, mf_g);
    // Compute L2 norm of the error
    const lf::fe::MeshFunctionFE mf_uh(fes_p, mu_vec);
    const double L2_err = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_uex - mf_uh), 8));
    errs.emplace_back(mu_vec.size(), L2_err);
  }
  // Output table of errors to file and terminal
  std::cout << "FE space type = " << typeid(LAGRFESPACE).name() << std::endl;
  std::ofstream out_file("errors.txt");
  std::cout.precision(3);
  std::cout << std::endl
            << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 err\n";
  std::cout << "---------------------------------------------" << '\n';
  for (const auto& err : errs) {
    auto [N, L2err] = err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << L2err << '\n';
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << L2err << '\n';
  }
}

}  // namespace LeastSquaresAdvection

#endif
