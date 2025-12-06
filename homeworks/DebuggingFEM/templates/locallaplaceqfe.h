/**
 * @file locallaplaceqfe.h
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NPDECODES_DEBUGGINGFEM_LOCALLAPLACEQFE_H_
#define NPDECODES_DEBUGGINGFEM_LOCALLAPLACEQFE_H_

#include <lf/mesh/mesh.h>

#include <Eigen/Core>

namespace DebuggingFEM {

class EntityMatrixProvider {
 public:
  virtual ~EntityMatrixProvider() = default;
  virtual bool isActive(const lf::mesh::Entity& cell) = 0;
  virtual Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity& cell) = 0;
};

class LocalLaplaceQFE1 : public EntityMatrixProvider {
 public:
  bool isActive(const lf::mesh::Entity& cell) override { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity& cell) override;
};

class LocalLaplaceQFE2 : public EntityMatrixProvider {
 public:
  bool isActive(const lf::mesh::Entity& cell) override { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity& cell) override;
};

class LocalLaplaceQFE3 : public EntityMatrixProvider {
 public:
  bool isActive(const lf::mesh::Entity& cell) override { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity& cell) override;
};

}  // namespace DebuggingFEM

#endif
