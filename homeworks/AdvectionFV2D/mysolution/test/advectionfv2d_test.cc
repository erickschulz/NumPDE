/**
 * @file advectionfv2d_test.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 13.07.2020
 * @copyright Developed at ETH Zurich
 */

#include "../advectionfv2d.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace AdvectionFV2D::test {

TEST(AdvectionFV2D, computeCellNormals) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);

  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      nv = AdvectionFV2D::computeCellNormals(mesh_p);

  Eigen::MatrixXd elem0(2, 3);
  elem0 << 0.55470, -1.0, 0.5547, 0.83205, 0.0, -0.83205;
  Eigen::MatrixXd elem1(2, 3);
  elem1 << -0.0, 0.894427, -0.707107, -1.0, 0.447214, 0.707107;
  Eigen::MatrixXd elem2(2, 3);
  elem2 << 0.894427, -0.0, -0.894427, -0.447214, 1.0, -0.447214;
  Eigen::MatrixXd elem3(2, 3);
  elem3 << -0.0, 0.707107, -0.894427, -1.0, 0.707107, 0.447214;
  Eigen::MatrixXd elem4(2, 3);
  elem4 << 1.0, -0.707107, -0.707107, -0.0, 0.707107, -0.707107;
  Eigen::MatrixXd elem5(2, 4);
  elem5 << 0.707107, 0.894427, -0.5547, -1.0, -0.707107, -0.447214, 0.83205,
      -0.0;
  Eigen::MatrixXd elem6(2, 4);
  elem6 << 0.707107, 1.0, -0.5547, -0.894427, -0.707107, -0.0, 0.83205,
      -0.447214;
  Eigen::MatrixXd elem7(2, 3);
  elem7 << 0.5547, -0.0, -0.5547, -0.83205, 1.0, -0.83205;
  Eigen::MatrixXd elem8(2, 3);
  elem8 << -0.894427, -0.0, 0.894427, 0.447214, -1.0, 0.447214;

  double tol = 1.0e-6;
  auto el = mesh_p->Entities(0);

  ASSERT_NE(nv, nullptr);
  ASSERT_NEAR(0.0, (elem0 - (*nv)(*el[0])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem1 - (*nv)(*el[1])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem2 - (*nv)(*el[2])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem3 - (*nv)(*el[3])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem4 - (*nv)(*el[4])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem5 - (*nv)(*el[5])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem6 - (*nv)(*el[6])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem7 - (*nv)(*el[7])).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (elem8 - (*nv)(*el[8])).lpNorm<Eigen::Infinity>(), tol);
}

TEST(AdvectionFV2D, getAdjacentCellPointers) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);

  std::shared_ptr<
      lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
      ac = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  auto el = mesh_p->Entities(0);

  std::array<const lf::mesh::Entity*, 4> elem0{
      {el[7], nullptr, el[5], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem1{
      {nullptr, el[2], el[5], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem2{{el[3], el[8], el[1], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem3{
      {nullptr, el[4], el[2], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem4{
      {nullptr, el[6], el[3], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem5{{el[1], el[8], el[0], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem6{{el[4], nullptr, el[7], el[8]}};
  std::array<const lf::mesh::Entity*, 4> elem7{
      {el[6], nullptr, el[0], nullptr}};
  std::array<const lf::mesh::Entity*, 4> elem8{{el[5], el[2], el[6], nullptr}};

  ASSERT_NE(ac, nullptr);
  ASSERT_EQ(elem0, (*ac)(*el[0]));
  ASSERT_EQ(elem1, (*ac)(*el[1]));
  ASSERT_EQ(elem2, (*ac)(*el[2]));
  ASSERT_EQ(elem3, (*ac)(*el[3]));
  ASSERT_EQ(elem4, (*ac)(*el[4]));
  ASSERT_EQ(elem5, (*ac)(*el[5]));
  ASSERT_EQ(elem6, (*ac)(*el[6]));
  ASSERT_EQ(elem7, (*ac)(*el[7]));
  ASSERT_EQ(elem8, (*ac)(*el[8]));
}

////////////////////////////////////////////////////////////////////////////////
// TODO: Here I ignored the boundary condition (general case)
// Note in problem description or adapt this test
////////////////////////////////////////////////////////////////////////////////
TEST(AdvectionFV2D, initializeMOLODEMatrix) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);

  const lf::assemble::UniformFEDofHandler dofh(
      mesh_p, {{lf::base::RefEl::kPoint(), 0},
               {lf::base::RefEl::kSegment(), 0},
               {lf::base::RefEl::kTria(), 1},
               {lf::base::RefEl::kQuad(), 1}});

  Eigen::MatrixXd B(dofh.NumDofs(), dofh.NumDofs());
  B << -1.88562, 0, 0, 0, 0, 1.23744, 0, 0.648181, 0, 0, -1.06066, 0, 0, 0, 0,
      0, 0, 0, 0, 0.176777, -2.12132, 1.94454, 0, 0, 0, 0, 0, 0, 0, 0, -3.18198,
      0, 0, 0, 0, 0, 0, 0, 0, 1.41421, -2.82843, 0, 0, 0, 0, 0, 0.707107, 0, 0,
      0, -2.20971, 0, 0, 1.5026, 0, 0, 0, 0, 1.88562, 0, -3.06413, 0, 0, 0, 0,
      0, 0, 0, 0, 2.7695, -2.7695, 0, 0, 0, 2.12132, 0, 0, 0, 0.883883, 0,
      -3.0052;

  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      normal_vectors = AdvectionFV2D::computeCellNormals(mesh_p);

  std::shared_ptr<
      lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  Eigen::MatrixXd B_matrix = AdvectionFV2D::initializeMOLODEMatrix(
      dofh, beta, adjacentCells, normal_vectors);

  double tol = 1.0e-5;
  ASSERT_NEAR(0.0, (B_matrix - B).lpNorm<Eigen::Infinity>(), tol);
}

TEST(AdvectionFV2D, computeHmin) {
  std::array<double, 4> hmin_ref{{0.222222, 0.100154, 0.0500771, 0.0250386}};

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);

  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 3)};

  double tol = 1.0e-6;
  int num_meshes = mesh_seq_p->NumLevels();
  for (int level = 0; level < num_meshes; ++level) {
    auto cur_mesh = mesh_seq_p->getMesh(level);

    double hmin = AdvectionFV2D::computeHmin(cur_mesh);
    ASSERT_NEAR(hmin_ref.at(level), hmin, tol);
  }
}

TEST(AdvectionFV2D, simulateAdvection) {
  double T = 1.0;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);
  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 2)};

  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  Eigen::Vector2d x0(0.8, 0.2);
  double d = 0.2;
  auto u0 = [x0, d](Eigen::Vector2d x) -> double {
    double dist = (x - x0).norm();
    if (dist < d) {
      return std::pow(std::cos(M_PI / (2.0 * d) * dist), 2);
    }
    return 0.0;
  };

  auto cur_mesh = mesh_seq_p->getMesh(2);
  const lf::assemble::UniformFEDofHandler cur_dofh(
      cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                 {lf::base::RefEl::kSegment(), 0},
                 {lf::base::RefEl::kTria(), 1},
                 {lf::base::RefEl::kQuad(), 1}});

  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      normal_vectors = AdvectionFV2D::computeCellNormals(cur_dofh.Mesh());

  std::shared_ptr<
      lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(cur_dofh.Mesh());

  Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(
      cur_dofh, beta, u0, adjacentCells, normal_vectors, T);

  Eigen::VectorXd ref_res(cur_dofh.NumDofs());
  ref_res << 0.121615, 0.0834347, 0.0613004, 0.0987587, 0, 0.0447353, 0.0129572,
      0.0289471, 0, 0.0226698, 0.00735303, 0.0122238, 0.0590068, 0.0171464,
      0.041362, 0.0296812, 0, 0, 0, 0, 0, 0, 0, 0, 6.05675e-05, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0.0110078, 0.00120791, 0.00749582, 0.00337274,
      0.000311443, 0, 0.00131884, 0.000267453, 0, 0.000839601, 0, 0, 0, 0, 0, 0,
      0, 0, 4.57313e-06, 3.60971e-07, 0.00390858, 0.000346692, 0.000315199,
      0.00128013, 1.32502e-06, 6.84524e-05, 7.24485e-05, 1.60551e-05, 0, 0,
      0.000279443, 2.38696e-05, 0, 0, 0.00365695, 0, 0.0291664, 0.00476745,
      0.0204668, 0.0155658, 3.43396e-08, 0.00641576, 0.00179795, 0.00254441, 0,
      0.00039036, 0.00427955, 0.00105022, 0.00262051, 0.0243687, 0.0113867,
      0.00115973, 0.130722, 0.0806817, 0.0410473, 0.0768465, 0.0105664,
      0.0354685, 0.016968, 0, 0.0744711, 0.0339141, 0.0450382, 0.119283, 0, 0,
      0.00504242, 0.0045039, 0, 0, 0.00520298, 0.00499163, 0.151107, 0.051108,
      0.0510654, 0.147083, 0.164793, 0.0643043, 0.142986, 0.161815, 0,
      0.0109585, 0.0123859, 0.0119805, 0.0250966, 0.0748938, 0.0498637,
      0.0424264, 0.0702418, 0.0602004, 0.104165, 0.0703989, 0.185889, 0.142512,
      0.1603, 0.184976, 0.00527327, 0.0475363, 0.01198, 0.0144631, 0.0300208,
      0.0936578, 0.0232557, 0.0515832, 0.0923879, 0.0344627, 0.125253,
      0.0774075;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (result - ref_res).lpNorm<Eigen::Infinity>(), tol);
}

TEST(AdvectionFV2D, findCFLthreshold) {
  double T = 1.0;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);
  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 4)};

  auto cur_mesh = mesh_seq_p->getMesh(4);

  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  const lf::assemble::UniformFEDofHandler cur_dofh(
      cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                 {lf::base::RefEl::kSegment(), 0},
                 {lf::base::RefEl::kTria(), 1},
                 {lf::base::RefEl::kQuad(), 1}});

  int threshold = AdvectionFV2D::findCFLthreshold(cur_dofh, beta, T);

  double tol = 3;
  ASSERT_NEAR(53, threshold, tol);
}

}  // namespace AdvectionFV2D::test
