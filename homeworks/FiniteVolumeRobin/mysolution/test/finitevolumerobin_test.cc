/**
 * @file finitevolumerobin_test.cc
 * @brief NPDE homework FiniteVolumeRobin code
 * @author Philippe Peter
 * @date February 2020
 * @copyright Developed at ETH Zurich
 */

#include "../finitevolumerobin.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

namespace FiniteVolumeRobin::test {

TEST(FiniteVolumeRobin, EdgeMatrixProvider) {
  // Build triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // initialize functors gamma and v
  auto gamma = [](Eigen::Vector2d x) { return 2.0; };
  auto v = [](Eigen::Vector2d x) { return 2 * x(0) + x(1); };
  auto v_mf = lf::mesh::utils::MeshFunctionGlobal(v);

  // set up finite element space and dofhandler
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto& dofh = fe_space->LocGlobMap();

  // mark boudnary edges
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // assemble galerkin matrix corresponding to correction term.
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  EdgeMatrixProvider edmat_provider(gamma, bd_flags);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edmat_provider, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // project v into the FEspace
  auto v_vec = lf::fe::NodalProjection<double>(*fe_space, v_mf);

  // create a ones vector
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(dofh.NumDofs());

  // evaluated product should corresponds to the integral of v*gamma over the
  // boundary of the mesh. Explanation: \union_i(\partial C_i \cap \partial
  // \Omega) = \partial \Omega, and the union is disjoint.
  auto product = (ones.transpose() * A_crs * v_vec).eval();
  EXPECT_NEAR(product(0, 0), 108, 1E-6);
}

TEST(FiniteVolumeRobin, EdgeVectorProvider) {
  // Build triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // initialize functor g
  auto g = [](Eigen::Vector2d x) { return 2 * x(0) + x(1); };

  // set up finite element space and dofhandler
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto& dofh = fe_space->LocGlobMap();

  // mark boudnary edges
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // assemble rhs vector
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();
  EdgeVectorProvider edvec_provider(g, bd_flags);
  lf::assemble::AssembleVectorLocally(1, dofh, edvec_provider, phi);

  // create a ones vector
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(dofh.NumDofs());

  // evaluated product should corresponds to the integral of g over the boundary
  // of the mesh. Explanation: \union_i(\partial C_i \cap \partial \Omega) =
  // \partial \Omega, and the union is disjoint
  auto product = (ones.transpose() * phi).eval();
  EXPECT_NEAR(product(0, 0), 54.0, 1E-6);
}

}  // namespace FiniteVolumeRobin::test
