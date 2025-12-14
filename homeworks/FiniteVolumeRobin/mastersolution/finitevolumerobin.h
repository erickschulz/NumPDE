/**
 * @file finitevolumerobin.h
 * @brief NPDE homework FiniteVolumeRobin code
 * @author Philippe Peter
 * @date February 2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>

namespace FiniteVolumeRobin {

/**
 * @headerfile finitevolumerobin.h
 * @brief Computes the local matrices on an edge which modify the "standard"
 * Galerkin matrix in the finite volume scheme in the exercise. Computations are
 * based on linear finite elements and the trapezoidal rule.
 *
 * @tparam FUNCTOR function that defines the scalar valued coefficient  \f$
 * \gamma \f$
 * @tparam EDGESELECTOR predicate defining which edges are included
 *
 * This helper class corresponds to the local matrix for the form
 * @f[
 *     \int\limits_{e \cap \partial C_i} \gamma(x) b^j(x),\mathrm{d}S(x)\;,
 * @f]
 *
 * where $@fe$@f is a (boundary) edge of the mesh and @f$ gamma \f$ is a
 * scalar-valued coefficient function.
 */
template <typename FUNCTOR, typename EDGESELECTOR>
class EdgeMatrixProvider {
 public:
  /**
   * @brief
   * @param gamma coefficient function
   * @param edge_selector predicate object selecting active edges to be covered
   * in the assembly (for the finite volume scheme considered in the exercise
   * only the boundary edges of the mesh are active)
   */
  explicit EdgeMatrixProvider(FUNCTOR gamma, EDGESELECTOR edge_selector)
      : gamma_(gamma), edge_sel_(edge_selector) {}

  /**
   * @brief actual computation of the edge matrix
   * @param edge reference to the edge for which the matrix is needed
   * @return a 2x2 dense matrix, contatining the edge matrix.
   *
   * Actual computation is based on linear finite elements and the trapezoidal
   * rule.
   */
  Eigen::Matrix2d Eval(const lf::mesh::Entity& edge);

  /**
   * @brief If true, then an edge is taken into account during assembly
   *
   * The information about "active" edges is supplied through the
   * `edge_selector` argument of the constructor.
   *
   * (For the finite volume scheme considered in the exercise only the boundary
   * edges are active.)
   */
  bool isActive(const lf::mesh::Entity& edge) const {
    LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for an edge");
    return edge_sel_(edge);
  }

 private:
  FUNCTOR gamma_;          // functor for the coefficient
  EDGESELECTOR edge_sel_;  // defines the active edges
};

template <typename FUNCTOR, typename EDGESELECTOR>
Eigen::Matrix2d EdgeMatrixProvider<FUNCTOR, EDGESELECTOR>::Eval(
    const lf::mesh::Entity& edge) {
  LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                "Function only defined on segments");
  const lf::geometry::Geometry* geo_ptr = edge.Geometry();
  Eigen::Matrix2d loc_mat;

  // get the length of the edge
  const double length = lf::geometry::Volume(*geo_ptr);

  // calculate the midpoint of the edge
  const Eigen::Matrix2d endpoints = lf::geometry::Corners(*geo_ptr);
  const Eigen::Vector2d midpoint = 0.5 * (endpoints.col(0) + endpoints.col(1));

  // compute the local matrix using the trapezoidal rule
  loc_mat(0, 0) =
      length * (gamma_(endpoints.col(0)) / 4.0 + gamma_(midpoint) / 8.0);
  loc_mat(1, 1) =
      length * (gamma_(endpoints.col(1)) / 4.0 + gamma_(midpoint) / 8.0);
  loc_mat(0, 1) = length * gamma_(midpoint) / 8.0;
  loc_mat(1, 0) = length * gamma_(midpoint) / 8.0;
  return loc_mat;
}

/**
 * @headerfile finitevolumerobin.h
 * @brief Computes the local vectors on an edge for the rhs of
 * the finite volume scheme in the exercise. Computations are based on the
 * trapezoidal rule.
 *
 * @tparam FUNCTOR function that defines the scalar valued coefficient  \f$ g
 * \f$
 * @tparam EDGESELECTOR predicate defining which edges are included
 *
 * This helper class corresponds to the local vector for the form
 * @f[
 *     \int\limits_{e \cap \partial C_i} g(x) ,\mathrm{d}S(x)\;,
 * @f]
 *
 * where $@fe$@f is a (boundary) edge of the mesh and @f$ g \f$ is a
 * scalar-valued coefficient function.
 */
template <typename FUNCTOR, typename EDGESELECTOR>
class EdgeVectorProvider {
 public:
  /**
   * @brief
   * @param g coefficient function
   * @param edge_selector predicate object selecting active edges to be covered
   * in the assembly (for the finite volume scheme considered in the exercise
   * only the boundary edges of the mesh are active)
   */
  explicit EdgeVectorProvider(FUNCTOR g, EDGESELECTOR edge_selector)
      : g_(g), edge_sel_(edge_selector) {}

  /**
   * @brief actual computation of the edge vector
   * @param edge reference to the edge for which the vector is needed
   * @return a 2x1 dense vector, contatining the edge matrix.
   *
   * Actual computation is based on the trapezoidal rule.
   */
  Eigen::Vector2d Eval(const lf::mesh::Entity& edge);

  /**
   * @brief If true, then an edge is taken into account during assembly
   *
   * The information about "active" edges is supplied through the
   * `edge_selector` argument of the constructor.
   *
   * (For the finite volume scheme considered in the exercise only the boundary
   * are active.)
   */
  bool isActive(const lf::mesh::Entity& edge) const {
    LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for an edge");
    return edge_sel_(edge);
  }

 private:
  FUNCTOR g_;              // functor for the coefficient
  EDGESELECTOR edge_sel_;  // defines the active edges
};

template <typename FUNCTOR, typename EDGESELECTOR>
Eigen::Vector2d EdgeVectorProvider<FUNCTOR, EDGESELECTOR>::Eval(
    const lf::mesh::Entity& edge) {
  LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                "Function only defined on segments");
  const lf::geometry::Geometry* geo_ptr = edge.Geometry();
  Eigen::Vector2d loc_vec;

  // get the length of the edge
  const double length = lf::geometry::Volume(*geo_ptr);

  // calculate the midpoint of the edge
  const Eigen::Matrix2d endpoints = lf::geometry::Corners(*geo_ptr);
  const Eigen::Vector2d midpoint = 0.5 * (endpoints.col(0) + endpoints.col(1));

  // compute the local vector using the trapezoidal rule
  loc_vec(0) = length / 4.0 * (g_(endpoints.col(0)) + g_(midpoint));
  loc_vec(1) = length / 4.0 * (g_(endpoints.col(1)) + g_(midpoint));
  return loc_vec;
}

}  // namespace FiniteVolumeRobin
