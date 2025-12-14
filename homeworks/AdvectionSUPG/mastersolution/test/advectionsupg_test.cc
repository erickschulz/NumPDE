/**
 * @file advectionsupg_test.cc
 * @brief NPDE homework AdvectionSUPG code
 * @author R. Hiptmair
 * @date July 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../advectionsupg.h"

#include <gtest/gtest.h>
#include <lf/assemble/assembler.h>
#include <lf/base/base.h>
#include <lf/base/ref_el.h>
#include <lf/fe/fe_tools.h>
#include <lf/geometry/geometry_interface.h>
#include <lf/mesh/utils/mesh_function_global.h>
#include <lf/quad/quad_rule.h>

#include <iostream>

#include "lf/mesh/test_utils/test_meshes.h"

namespace AdvectionSUPG::test {

// Output element matrices
template <class ELEMAT_PROVIDER>
void printElementMatrices(const lf::mesh::Mesh& mesh,
                          ELEMAT_PROVIDER& elmat_builder) {
  // Traverse the cells of the mesh and compute element matrices
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    const typename ELEMAT_PROVIDER::ElemMat M{elmat_builder.Eval(*cell)};
    const Eigen::MatrixXd vertices = lf::geometry::Corners(*cell->Geometry());
    std::cout << "=== Element matrix for triangle\n " << vertices << ":"
              << std::endl
              << M << std::endl;
  }
}

template <typename V_FUNCTOR, typename U_FUNCTOR, typename W_FUNCTOR,
          typename GRAD_U_FUNCTOR, typename GRAD_W_FUNCTOR>
bool checkBilinearForm(V_FUNCTOR&& v, U_FUNCTOR&& u, W_FUNCTOR&& w,
                       GRAD_U_FUNCTOR&& grad_u, GRAD_W_FUNCTOR&& grad_w) {
  // I. Direct integration
  // Build a lambda function for the integrand
  auto itg = [&](Eigen::Vector2d x) -> double {
    auto v_val = v(x);
    return v_val.dot(grad_u(x)) * (v_val.dot(grad_w(x)) + w(x));
  };
  // Perform numerical quadrature over the unit square using a high-order
  // Gaussian quadrature rule
  const lf::quad::QuadRule tpqr{
      lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 20)};
  const lf::base::size_type P = tpqr.NumPoints();
  const Eigen::MatrixXd zeta{tpqr.Points()};
  const Eigen::VectorXd omega{tpqr.Weights()};
  double val_itg = 0.0;
  for (unsigned int l = 0; l < P; ++l) {
    val_itg += omega[l] * itg(zeta.col(l));
  }

  // II. Inegration by means of bilinear form
  // Wrap functors into MeshFunctions
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  lf::mesh::utils::MeshFunctionGlobal mf_w(w);
  lf::mesh::utils::MeshFunctionGlobal mf_v(v);
  // Generate triangular mesh of the unit square
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Set up finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_p);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NumDofs());
  // Interpolate the functions u and w
  auto u_vec{lf::fe::NodalProjection(*fe_space, mf_u)};
  auto w_vec{lf::fe::NodalProjection(*fe_space, mf_w)};
  // Assemble Galerkin matrix
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  AdvectionSUPG::SUAdvectionElemMatrixProvider elmat_builder(mf_v, false);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  const double val_blf = w_vec.dot(A.MatVecMult(1.0, u_vec));
  std::cout << "val_itg = " << val_itg << " <-> val_blf = " << val_blf
            << std::endl;
  return ((std::abs(val_itg - val_blf) < 1.0E-6 * std::abs(val_itg)) ||
          (std::abs(val_itg - val_blf) < 1.0E-8));
}

TEST(AdvectionSUPG, printmat) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Velocity field
  auto vf = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(x[0], x[1]);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_v(vf);
  AdvectionSUPG::SUAdvectionElemMatrixProvider elmat_builder(mf_v);
  // Output element matrices
  printElementMatrices(*mesh_p, elmat_builder);
}

TEST(AdvectionSUPG, constfun) {
  // Test with constant functions
  auto u = [](Eigen::Vector2d /*x*/) -> double { return 1.0; };
  auto w = [](Eigen::Vector2d x) -> double { return x[0] + 2 * x[1]; };
  auto grad_u = [](Eigen::Vector2d /*x*/) -> Eigen::Vector2d {
    return Eigen::Vector2d(0.0, 0.0);
  };
  auto grad_w = [](Eigen::Vector2d /*x*/) -> Eigen::Vector2d {
    return Eigen::Vector2d(1.0, 2.0);
  };
  auto v = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  EXPECT_TRUE(checkBilinearForm(v, u, w, grad_u, grad_w));
}

TEST(AdvectionSUPG, linfun) {
  // Test with constant functions
  auto u = [](Eigen::Vector2d x) -> double { return (x[0] - x[1]); };
  auto w = [](Eigen::Vector2d x) -> double { return x[0] + 2 * x[1]; };
  auto grad_u = [](Eigen::Vector2d /*x*/) -> Eigen::Vector2d {
    return Eigen::Vector2d(1.0, -1.0);
  };
  auto grad_w = [](Eigen::Vector2d /*x*/) -> Eigen::Vector2d {
    return Eigen::Vector2d(1.0, 2.0);
  };
  auto v = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  EXPECT_TRUE(checkBilinearForm(v, u, w, grad_u, grad_w));
}

TEST(AdvectionSUPG, quadfun) {
  // Test with constant functions
  auto u = [](Eigen::Vector2d x) -> double { return (x[0] * x[1] - x[1]); };
  auto w = [](Eigen::Vector2d x) -> double { return x[0] + 2.0 * x[1]; };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(x[1], x[0] - 1.0);
  };
  auto grad_w = [](Eigen::Vector2d /*x*/) -> Eigen::Vector2d {
    return Eigen::Vector2d(1.0, 2.0);
  };
  auto v = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  EXPECT_TRUE(checkBilinearForm(v, u, w, grad_u, grad_w));
}

}  // namespace AdvectionSUPG::test
