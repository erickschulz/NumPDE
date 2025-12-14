/**
 * @file
 * @brief NPDE homework ResidualErrorEstimator
 * @author Ralf Hiptmair
 * @date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <iostream>
#include <memory>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace REE {

/** @brief MeshFunction type representing a real-valued piecewise constant
 * function
 *
 * This is rudimentary implementation tha merely allows initialization through a
 * function object
 *
 * @note Evaluation is possible only for cells!
 */
class MeshFunctionPWConst {
 public:
  /** @brief Constructor sampling a function at the centers of cells
   *
   * @tparam FUNCTOR a type providing an Eigen::Vector2d -> double evaluation
   * operator
   * @param mesh_p shared pointer to the mesh the MeshFunction belongs to
   * @param f the function object, with which the MeshFunction is initialized
   */
  template <typename FUNCTOR>
  MeshFunctionPWConst(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                      FUNCTOR&& f);
  /** @brief evaluation operator
   *
   * @param e reference to a cell of the underlying mesh
   * @param refc reference coordinates of evaluation points. Only their number
   * is used!
   * @return vector of function values, all the same!
   */
  std::vector<double> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& refc) const {
    LF_ASSERT_MSG(e.RefEl().Dimension() == 2, "Implemented for 2D cells only!");
    return std::vector<double>(refc.cols(), data_(e));
  }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  lf::mesh::utils::CodimMeshDataSet<double> data_;
};

template <typename FUNCTOR>
MeshFunctionPWConst::MeshFunctionPWConst(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, FUNCTOR&& f)
    : mesh_p_(mesh_p), data_(mesh_p, 0, 0.0) {
  // Run through the cells of the mesh
  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    // Fetch type of cell
    const lf::base::RefEl ref_el{cell->RefEl()};
    // Get geometry, corners of cell
    const lf::geometry::Geometry& geo{*(cell->Geometry())};
    const Eigen::MatrixXd corners{lf::geometry::Corners(geo)};
    // Compute physical coordinates of center point
    Eigen::MatrixXd center(2, 1);
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        center = corners * Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        break;
      }
      case lf::base::RefEl::kQuad(): {
        center = corners * Eigen::Vector4d(0.25, 0.25, 0.25, 0.25);
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Unknown cell type");
      }
    }
    // Evaluate given function. Return value must allow conversion to double!
    data_(*cell) = f(center);
  }
}

/** @brief Data describing the FE-discretized boundary-value problem

    Theses data comprise a linear Lagrangian finite-element space, the diffusion
    coefficient as a MeshFunction, and the right-hand side source function as
   another MeshFunction.
 */
/* SAM_LISTING_BEGIN_1 */
struct dataDiscreteBVP {
  dataDiscreteBVP& operator=(const dataDiscreteBVP&) = delete;
  dataDiscreteBVP& operator=(const dataDiscreteBVP&&) = delete;

  /** @brief Constructor, which essentially copies the passed arguments
   */
  dataDiscreteBVP(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                  std::function<double(Eigen::Vector2d)> alpha,
                  std::function<double(Eigen::Vector2d)> f);
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> pwlinfespace_p_;
  lf::mesh::utils::MeshFunctionGlobal<std::function<double(Eigen::Vector2d)>>
      mf_f_;
  MeshFunctionPWConst mf_alpha_;
};
/* SAM_LISTING_END_1 */

/** @briefs Solves homogeneous Dirichlet boundary value problem
 */
Eigen::VectorXd solveBVP(const dataDiscreteBVP& disc_bvp);

/** @brief Computes cell contributions to error estimator */
lf::mesh::utils::CodimMeshDataSet<double> volumeResiduals(
    const dataDiscreteBVP& disc_bvp, const Eigen::VectorXd& u_vec);

/** @brief Evaluates edge terms for  error estimator */
lf::mesh::utils::CodimMeshDataSet<double> edgeResiduals(
    const dataDiscreteBVP& disc_bvp, const Eigen::VectorXd& u_vec);

/** @brief solves boundary value problem and estimates error

    @note the provided mesh must cover the unit square
 */
std::tuple<double, double, double> solveAndEstimate(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

}  // namespace REE
