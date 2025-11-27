/**
 * @file taylorhoodnonmonolithic.cc
 * @brief NPDE homework TaylorHoodNonMonolithic code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "taylorhoodnonmonolithic.h"

#include <lf/base/lf_assert.h>
#include <lf/fe/fe_tools.h>
#include <lf/fe/mesh_function_fe.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/utils/mesh_function_unary.h>

namespace TaylorHoodNonMonolithic {
/* SAM_LISTING_BEGIN_1 */
THBElementMatrixProvider::ElemMat THBElementMatrixProvider::Eval(
    const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << cell.RefEl());
  // Element matrix to be filled
  ElemMat BK;
  // Area of the triangle
  double area = lf::geometry::Volume(*cell.Geometry());
  // Compute gradients of barycentric coordinate functions, see
  // \lref{cpp:gradbarycoords}.
  // Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(cell.Geometry()));
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();
  // This matrix contains $\cob{\grad \lambda_i}$ in its columns
  const auto G{X.inverse().block<2, 3>(1, 0)};
  // Dummy lambda functions for barycentric coordinates
  std::array<std::function<double(Eigen::Vector3d)>, 3> lambda{
      [](Eigen::Vector3d c) -> double { return c[0]; },
      [](Eigen::Vector3d c) -> double { return c[1]; },
      [](Eigen::Vector3d c) -> double { return c[2]; }};
  // Gradients of local shape functions of quadratic Lagrangian finite element
  // space as lambda functions, see \prbeqref{eq:quadlsf}
  std::array<std::function<Eigen::Vector2d(Eigen::Vector3d)>, 6> gradbK{
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return (4 * c[0] - 1) * G.col(0);
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return (4 * c[1] - 1) * G.col(1);
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return (4 * c[2] - 1) * G.col(2);
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return 4 * (c[0] * G.col(1) + c[1] * G.col(0));
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return 4 * (c[1] * G.col(2) + c[2] * G.col(1));
      },
      [&G](Eigen::Vector3d c) -> Eigen::Vector2d {
        return 4 * (c[2] * G.col(0) + c[0] * G.col(2));
      }};
  // Barycentric coordinates of the midpoints of the edges for
  // use with the 3-point edge midpoint quadrature rule \prbeqref{eq:MPR}
  const std::array<Eigen::Vector3d, 3> mp = {Eigen::Vector3d({0.5, 0.5, 0}),
                                             Eigen::Vector3d({0, 0.5, 0.5}),
                                             Eigen::Vector3d({0.5, 0, 0.5})};
#if SOLUTION
  // Initialize the element matrix for the bilinear form b(.,.)
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      // \prbeqref{eq:BKtrp} without the area scaling
      const Eigen::Vector2d gql_ij{gradbK[j](mp[0]) * lambda[i](mp[0]) +
                                   gradbK[j](mp[1]) * lambda[i](mp[1]) +
                                   gradbK[j](mp[2]) * lambda[i](mp[2])};
      BK(i, j) = gql_ij[dirflag_];
    }
  }
  // Do not forget to scale by $\frac13|K|$.
  BK *= area / 3.0;
#else
  /* **********************************************************************
     Your code here
     ********************************************************************** */
#endif
  return BK;
}
/* SAM_LISTING_END_1 */
#if SOLUTION
/** @brief Monitor Uzawa convergende
 *
 */
Eigen::Matrix<double, Eigen::Dynamic, 6> monitorUzawaConvergence(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, bool print) {
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!i
   * This function assumes that the local -> global index mappings set up when
   * creating a finite element space for a given mesh are always the same
   * provided that the mesh object is the same.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  // Set up finite element space
  auto fes_LO2 =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh_LO2{fes_LO2->LocGlobMap()};
  [[maybe_unused]] const lf::assemble::size_type n_LO2 = dofh_LO2.NumDofs();
  auto fes_LO1 =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh_LO1{fes_LO1->LocGlobMap()};
  [[maybe_unused]] const lf::assemble::size_type n_LO1 = dofh_LO1.NumDofs();

  // Force functor (rotational force)
  auto force = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return {-x[1], x[0]};
  };
  // Generate Galerkin matrices
  auto [A, B_x, B_y, phi_x, phi_y] = buildStokesLSE(mesh_p, force);
  // For recording progress of the iteration
  std::vector<std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>>
      rec_data;
  auto rec = [&rec_data](const Eigen::VectorXd& mu_x,
                         const Eigen::VectorXd& mu_y,
                         const Eigen::VectorXd& pi) -> void {
    rec_data.emplace_back(mu_x, mu_y, pi);
  };
  auto [vec_res_mu_x, vec_res_mu_y, vec_res_pi] =
      TaylorHoodNonMonolithic::CGUzawa(A, B_x, B_y, phi_x, phi_y, 1E-8, 1E-10,
                                       100, rec);
  if (print) {
    std::cout << "\t >>>>> CG Uzawa took " << rec_data.size() << "steps\n";
  }

  // Matrix for returning iteration error norms
  Eigen::Matrix<double, Eigen::Dynamic, 6> it_err_norms(rec_data.size(), 6);
  it_err_norms.setZero();

  // Set up mesh functions for solution returned by iteration
  const lf::fe::MeshFunctionFE mf_res_mu_x(fes_LO2, vec_res_mu_x);
  const lf::fe::MeshFunctionFE mf_res_mu_y(fes_LO2, vec_res_mu_y);
  const lf::fe::MeshFunctionFE mf_res_pi(fes_LO1, vec_res_pi);
  const lf::fe::MeshFunctionGradFE gmf_res_mu_x(fes_LO2, vec_res_mu_x);
  const lf::fe::MeshFunctionGradFE gmf_res_mu_y(fes_LO2, vec_res_mu_y);
  int step = 1;
  for (auto vecs : rec_data) {
    auto& [mu_x, mu_y, pi] = vecs;
    const Eigen::VectorXd res_x = phi_x - A * mu_x - B_x.transpose() * pi;
    const Eigen::VectorXd res_y = phi_y - A * mu_y - B_y.transpose() * pi;
    const Eigen::VectorXd res_pi = B_x * mu_x + B_y * mu_y;
    const lf::fe::MeshFunctionFE mf_mu_x(fes_LO2, mu_x);
    const lf::fe::MeshFunctionFE mf_mu_y(fes_LO2, mu_y);
    const lf::fe::MeshFunctionFE mf_pi(fes_LO1, pi);
    const lf::fe::MeshFunctionGradFE gmf_mu_x(fes_LO2, mu_x);
    const lf::fe::MeshFunctionGradFE gmf_mu_y(fes_LO2, mu_y);
    it_err_norms(step - 1, 0) = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_mu_x - mf_res_mu_x), 4));
    it_err_norms(step - 1, 1) = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_mu_y - mf_res_mu_y), 4));
    it_err_norms(step - 1, 2) = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(gmf_mu_x - gmf_res_mu_x), 4));
    it_err_norms(step - 1, 3) = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(gmf_mu_y - gmf_res_mu_y), 4));
    it_err_norms(step - 1, 4) = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_pi - mf_res_pi), 4));
    it_err_norms(step - 1, 5) = res_pi.norm();
    if (print) {
      std::cout << "step " << step << ": |res_x| = " << res_x.norm()
                << ", |res_y| = " << res_y.norm()
                << ", |res_pi| = " << res_pi.norm() << std::endl;
      std::cout << "L2(mu_x) = " << it_err_norms(step - 1, 0);
      std::cout << ", L2(mu_y) = " << it_err_norms(step - 1, 1);
      std::cout << ", H1(mu_x) = " << it_err_norms(step - 1, 2);
      std::cout << ", H1(mu_y) = " << it_err_norms(step - 1, 3);
      std::cout << ", L2(pi) = " << it_err_norms(step - 1, 4) << std::endl;
    }
    step++;
  }
  return it_err_norms;
}
#endif

}  // namespace TaylorHoodNonMonolithic
