#ifndef POINTEVALUTION_H
#define POINTEVALUTION_H
/**
 * @ file pointEvaluation.h
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch, Liaowang Huang (refactoring)
 * @ date 22/03/2019, 06/01/2020 (refactoring)
 * @ copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <utility>

namespace PointEvaluationRhs {

std::pair<double, double> normsSolutionPointLoadDirichletBVP(
    const lf::assemble::DofHandler& dofh, Eigen::Vector2d source_point,
    Eigen::VectorXd& sol_vec);

Eigen::Vector2d GlobalInverseTria(Eigen::Matrix<double, 2, 3> mycorners,
                                  Eigen::Vector2d x);

Eigen::Vector2d GlobalInverseQuad(Eigen::Matrix<double, 2, 4> vert,
                                  Eigen::Vector2d x);

inline double triaArea(const Eigen::Vector2d a, const Eigen::Vector2d b,
                       const Eigen::Vector2d c);

std::pair<double, double> solveQuadraticEquation(double a, double b, double c);

class DeltaLocalVectorAssembler {
 private:
  Eigen::Vector2d x_0;
  bool already_found;

 public:
  explicit DeltaLocalVectorAssembler(Eigen::Vector2d x)
      : x_0(x), already_found(false) {}
  bool isActive(const lf::mesh::Entity& entity) const {
    return (!already_found);
  }
  Eigen::VectorXd Eval(const lf::mesh::Entity& cell);
};

}  // namespace PointEvaluationRhs
#endif  // define POINTEVALUATION_H
