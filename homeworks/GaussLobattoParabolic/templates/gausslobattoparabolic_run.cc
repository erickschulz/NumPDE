/**
 * @file gausslobattoparabolic_run.cc
 * @brief Main file for conducting convergence studies
 * @author Oliver Rietmann
 * @date 22.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "gausslobattoparabolic.h"

constexpr double PI = 3.14159265358979323846;

// Helper function to compute the meshwidth of a refinement
double maxLength(const std::span<const lf::mesh::Entity* const>& edges) {
  double length = 0.0;
  for (const lf::mesh::Entity* edge : edges) {
    Eigen::Matrix2d corners = lf::geometry::Corners(*(edge->Geometry()));
    double new_length = (corners.col(1) - corners.col(0)).norm();
    length = std::max(length, new_length);
  }
  return length;
}

int main() {
  // Load the mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  "meshes/hex_hybrid.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Solve the PDE
  double T = 1.0;
  unsigned int M = 100;
  auto g = [](double t) { return t < 1.0 ? std::sin(0.5 * PI * t) : 1.0; };
  auto u0 = [](Eigen::Vector2d x) { return 0.0; };
  Eigen::VectorXd mu =
      GaussLobattoParabolic::evolveIBVPGaussLobatto_u0(fe_space, T, M, g, u0);

  // Write the solution to a .vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "solution.vtk");
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  const lf::assemble::DofHandler& dofh = fe_space->LocGlobMap();
  for (int i = 0; i < mu.size(); ++i) {
    nodal_data->operator()(dofh.Entity(i)) = mu(i);
  };
  vtk_writer.WritePointData("solution", *nodal_data);

  // Compute specio-temporal convergence table

  // Build mesh hierarchy to compute the spatio-temporal discretization error
  std::shared_ptr<lf::mesh::Mesh> base_mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> mesh_hierarchy_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(base_mesh_p, 6);

  // Boundary data
  auto g_ = [](double t) { return 1.0; };

  // Inital data
  auto u0_ = [](Eigen::Vector2d x) {
    return std::sin(2.0 * PI * x(0)) * std::sin(2.0 * PI * x(1)) + 1.0;
  };

  /*// Exact solution at time t
  auto u = [](Eigen::Vector2d x, double t) {
    return std::exp((-8.0) * PI * PI * t) * std::sin(2.0 * PI * x(0)) *
  std::sin(2.0 * PI * x(1)) + 1.0;
  };*/

  // Gradient of the exact solution at time t
  auto grad_u = [](Eigen::Vector2d x, double t) -> Eigen::Vector2d {
    double amplitude = 2.0 * PI * std::exp((-8.0) * PI * PI * t);
    double x_ = 2.0 * PI * x(0);
    double y_ = 2.0 * PI * x(1);
    return amplitude * Eigen::Vector2d(std::cos(x_) * std::sin(y_),
                                       std::sin(x_) * std::cos(y_));
  };

  // Compute meshwidth for each refinement
  int L = mesh_hierarchy_p->NumLevels();
  Eigen::VectorXd meshwidth(L);
  Eigen::VectorXd spacetime_error(L);

  for (int k = 0; k < L; ++k) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = mesh_hierarchy_p->getMesh(k);
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Vectors for storing the H1-seminorm error at every timestep
    std::vector<double> h1s_error;

    // Fills the error vector that is captured by reference
    auto recorder = [&h1s_error, &fe_space, grad_u](double t,
                                                    const Eigen::VectorXd& mu) {
      lf::mesh::utils::MeshFunctionGlobal grad_mf(
          [t, grad_u](Eigen::Vector2d x) { return grad_u(x, t); });
      lf::uscalfe::MeshFunctionL2GradientDifference loc_comp_h1(fe_space,
                                                                grad_mf, 4);
      double error = lf::uscalfe::NormOfDifference(
          fe_space->LocGlobMap(), loc_comp_h1, mu, lf::base::PredicateTrue());
      h1s_error.push_back(error);
    };

    // Compute the meshwidth of the current refinement
    meshwidth(k) = maxLength(mesh_p->Entities(1));

    // Choose M such that the uniform timestep size \Blue{$\tau$}
    // is proportional to the square root of the meshwidth
    double T = 0.1;
    int M = (int)(10 * T * std::sqrt(1.0 / meshwidth(k)));
    GaussLobattoParabolic::evolveIBVPGaussLobatto(fe_space, T, M, g_, u0_,
                                                  recorder);

    // Sum the time-weighted H1-seminorm erros to get the spatio-temporal error
    Eigen::VectorXd error_at_t = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        h1s_error.data(), h1s_error.size());
    double tau = T / M;
    spacetime_error(k) = (std::sqrt(tau) * error_at_t).norm();
  }

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::cout << "meshwidth: " << meshwidth.transpose().format(CSVFormat)
            << std::endl;
  std::cout << "error:     " << spacetime_error.transpose().format(CSVFormat)
            << std::endl;
  return 0;
}
