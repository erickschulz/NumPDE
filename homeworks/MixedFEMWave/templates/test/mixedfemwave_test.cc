/**
 * @file mixedfemwave_test.cc
 * @author Erick Schulz
 * @date 24.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../mixedfemwave.h"

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

namespace MixedFEMWave::test {

TEST(MixedFEMWave_computeMQ, test) {
  // LOADING COARSE MESH
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 "meshes/unitsquare_unitest.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  // Vector dofhandler for the finite element space Q
  lf::assemble::UniformFEDofHandler dofh_Q(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 0},
                                            {lf::base::RefEl::kSegment(), 0},
                                            {lf::base::RefEl::kTria(), 2},
                                            {lf::base::RefEl::kQuad(), 2}});

  Eigen::SparseMatrix<double> M_Q_sparse = computeMQ(dofh_Q);
  Eigen::MatrixXd M_Q_dense = Eigen::MatrixXd(M_Q_sparse);
  Eigen::MatrixXd M_Q_dense_ref = Eigen::MatrixXd::Zero(8, 8);
  for (int i = 0; i < 8; i++) {
    M_Q_dense_ref(i, i) = 0.25;
  }

  double tol = 1.0e-8;
  Eigen::MatrixXd difference = M_Q_dense - M_Q_dense_ref;
  ASSERT_NEAR(0.0, difference.lpNorm<Eigen::Infinity>(), tol);
}

TEST(MixedFEMWave_computeMV, test) {
  auto rho = [](Eigen::Vector2d x) -> double { return 1.0; };
  // LOADING COARSE MESH
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 "meshes/unitsquare_unitest.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  // Scalar finite element space for lowest-order Lagrangian finite elements
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  Eigen::SparseMatrix<double> M_V_sparse = computeMV(fe_space_V, rho);
  Eigen::MatrixXd M_V_dense = Eigen::MatrixXd(M_V_sparse);
  Eigen::MatrixXd M_V_dense_ref = Eigen::MatrixXd::Zero(5, 5);
  for (int i = 0; i < 4; i++) {
    M_V_dense_ref(i, i) = 0.1 + 0.2 / 3.0;
  }
  M_V_dense_ref(4, 4) = 1.0 / 3.0;

  double tol = 1.0e-8;
  Eigen::MatrixXd difference = M_V_dense - M_V_dense_ref;
  ASSERT_NEAR(0.0, difference.lpNorm<Eigen::Infinity>(), tol);
}

TEST(MixedFEMWave_computeB, test) {
  auto rho = [](Eigen::Vector2d x) -> double { return 1.0; };
  // LOADING COARSE MESH
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 "meshes/unitsquare_unitest.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  // Scalar finite element space for lowest-order Lagrangian finite elements
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_V =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Scalar dofhandler as built along with the finite-element space for V
  const lf::assemble::DofHandler& dofh_V = fe_space_V->LocGlobMap();

  // Vector dofhandler for the finite element space Q
  lf::assemble::UniformFEDofHandler dofh_Q(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 0},
                                            {lf::base::RefEl::kSegment(), 0},
                                            {lf::base::RefEl::kTria(), 2},
                                            {lf::base::RefEl::kQuad(), 2}});

  Eigen::SparseMatrix<double> B_sparse = computeB(dofh_V, dofh_Q);
  Eigen::MatrixXd B_dense = Eigen::MatrixXd(B_sparse);
  Eigen::MatrixXd B_dense_ref = Eigen::MatrixXd::Zero(8, 5);
  B_dense_ref(1, 4) = 0.5;
  B_dense_ref(2, 4) = 0.5;
  B_dense_ref(4, 4) = -0.5;
  B_dense_ref(7, 4) = -0.5;

  double tol = 1.0e-8;
  Eigen::MatrixXd difference = B_dense - B_dense_ref;
  ASSERT_NEAR(0.0, difference.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace MixedFEMWave::test
