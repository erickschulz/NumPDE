/**
 * @file coupledbvps_test.cc
 * @brief NPDE homework CoupledBVPs code
 * @author W. Tonnon
 * @date June 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../coupledbvps.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <memory>

namespace CoupledBVPs::test {
TEST(CoupledBVPs, BVPsolver) {
  // Define the file that contains the mesh
  std::string mesh_file = "meshes/hexagon.msh";

  // Read and process the mesh file
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(factory), mesh_file);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  // Define a suitable mesh-hierarchy for convergence analysis
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 3);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  // We define an appropriate manufactured solution u with boundary condition g,
  // source term f, and coefficient alpha
  auto u_exact = [](Eigen::Vector2d x) -> double {
    return std::sin(M_PI * x(0)) * std::cos(M_PI * x(1));
  };
  auto g = [](Eigen::Vector2d x) -> double {
    return std::sin(M_PI * x(0)) * std::cos(M_PI * x(1));
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return M_PI * std::sin(M_PI * x(0)) * std::cos(M_PI * x(1));
  };
  auto alpha = [](Eigen::Vector2d x) -> double { return 1. / (2. * M_PI); };

  // We cast the lambda functions into meshfunction as needed by the
  // solveDirichletBVP function and error computation
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};
  lf::mesh::utils::MeshFunctionGlobal mf_u{u_exact};

  // We loop over the mesh hierarchy and check if the convergence rates are as
  // expected
  double L2err = 1.;
  for (int l = 0; l < 3; ++l) {
    // We extract the appropriate mesh and define the FE space on that mesh
    std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p = multi_mesh_p->getMesh(l);
    auto fes_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(lev_mesh_p);

    // We call solveDirichletBVP on the specified FE space
    Eigen::VectorXd sol =
        CoupledBVPs::solveDirichletBVP(fes_p, mf_alpha, mf_f, g);

    // We cast the solution obtained by solveDirichletBVP into a MeshFunction
    // for error computation
    const lf::fe::MeshFunctionFE mf_sol(fes_p, sol);

    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err_new = std::sqrt(lf::fe::IntegrateMeshFunction(
        *lev_mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));

    // Compare the error to the error on the previous mesh. For the L2 norm, we
    // expect second-order convergence
    EXPECT_GE(-std::log(L2err_new / L2err) / std::log(2.), 1.95);

    // Store the L2 error for the next iteration
    L2err = L2err_new;
  }
}

TEST(CoupledBVPs, solveModulatedHeatFlow2) {
  // Define the file that contains the mesh
  std::string mesh_file = "meshes/hexagon.msh";

  // Read and process the mesh file
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(factory), mesh_file);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  // Define a suitable mesh-hierarchy for convergence analysis
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 3);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  // We define a source term for w
  auto f = [](Eigen::Vector2d x) -> double { return std::sin(M_PI * x(0)); };

  // We define an FE space on the chosen mesh
  auto fes_p = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // We define a precomputed solution for w
  Eigen::VectorXd exact_w_gf(fes_p->LocGlobMap().NumDofs());
  exact_w_gf << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.03195, 1.03217,
      1.03464, 1.03922, 1.02975, 1.03052, 1.02957, 1.02928, 1.03712, 1.03803,
      1.02234, 1.02081, 1.03122, 1.03238, 1.05495, 1.05503, 1.07147, 1.07452,
      1.08347, 1.0842, 1.08887, 1.091, 1.09258, 1.09112, 1.09489, 1.08893,
      1.08453, 1.09316, 1.08967, 1.07614, 1.07297, 1.0853, 1.07598, 1.08822,
      1.08069, 1.06414, 1.0804, 1.06986, 1.08263, 1.06153, 1.07791, 1.07589,
      1.08296, 1.09309, 1.08009, 1.0891, 1.05901, 1.06908, 1.07614, 1.07662,
      1.0894, 1.0612, 1.06142, 1.0736, 1.053, 1.04977, 1.0577, 1.0886, 1.03374,
      1.03721, 1.03278, 1.05041, 1.02943, 1.02832, 1.03631, 1.04668, 1.07445,
      1.02049, 1.05138, 1.02349, 1.04465, 1.04421, 1.09284, 1.07447, 1.05483,
      1.08928, 1.07465, 1.08771, 1.01444, 1.01199, 1.05638, 1.05571, 1.08917,
      1.09493, 1.06479, 1.04611, 1.06161, 1.03529, 1.04256, 1.028, 1.022,
      1.02145, 1.05098, 1.05495, 1.06666, 1.06528, 1.02704, 1.06005, 1.06016,
      1.03041;

  // We define a precomputed solution for v
  Eigen::VectorXd exact_v_gf(fes_p->LocGlobMap().NumDofs());
  exact_v_gf << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.37132e-17, 0, 0, 0, 0, 0, 0.0534223,
      0.0537455, 0.0562894, 0.0634262, 0.0791405, 0.0811659, 0.0512667,
      0.0507815, 0.0608571, 0.0623526, 0.0416182, 0.0388594, 0.0546221,
      0.0732333, 0.131114, 0.133791, 0.170486, 0.167179, 0.193562, 0.189583,
      0.204747, 0.201117, 0.193752, 0.201632, 0.184047, 0.184767, 0.195531,
      0.170242, 0.155131, 0.170407, 0.173686, 0.192292, 0.164258, 0.173673,
      0.154781, 0.133593, 0.1347, 0.116145, 0.143184, 0.131335, 0.14832,
      0.140807, 0.144316, 0.16949, 0.134077, 0.153354, 0.118178, 0.114801,
      0.127027, 0.127892, 0.160207, 0.100727, 0.101488, 0.124836, 0.0911715,
      0.0836828, 0.114077, 0.183192, 0.0579527, 0.0831596, 0.074611, 0.0846228,
      0.0502402, 0.0483317, 0.0793333, 0.0771558, 0.135611, 0.0383241,
      0.0869284, 0.0434405, 0.07242, 0.0718252, 0.194224, 0.159972, 0.0921017,
      0.161346, 0.127262, 0.172296, 0.0288841, 0.0243543, 0.134809, 0.135154,
      0.205391, 0.183519, 0.113241, 0.0792201, 0.10547, 0.0583825, 0.0740272,
      0.0466403, 0.0388276, 0.0378929, 0.0940449, 0.100959, 0.112766, 0.109008,
      0.0437223, 0.099065, 0.0992796, 0.0551074;

  // We cast the precomputed solutions for w and v into MeshFunctions for error
  // computation
  lf::fe::MeshFunctionFE mf_w(fes_p, exact_w_gf);
  lf::fe::MeshFunctionFE mf_v(fes_p, exact_v_gf);

  // We call solveDirichletBVP on the specified FE space
  auto sol = CoupledBVPs::solveModulatedHeatFlow(fes_p, f);

  // We cast the solutions obtained by solveModulatedHeatFlow into a
  // MeshFunction for error computation
  const lf::fe::MeshFunctionFE mf_sol_w(fes_p, sol.first);
  const lf::fe::MeshFunctionFE mf_sol_v(fes_p, sol.second);

  // compute errors with 3rd order quadrature rules, which is sufficient for
  // piecewise linear finite elements
  double L2err_w = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(mf_sol_w - mf_w), 2));
  double L2err_v = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(mf_sol_v - mf_v), 2));

  // We expect the L2 integral of the error to be around machine precision
  // (around 1e-14), so the square root of this error should be around 1e-7. To
  // be safe, we ask for 1e-5.
  EXPECT_LT(L2err_w, 1e-5);
  EXPECT_LT(L2err_v, 1e-5);
}

}  // namespace CoupledBVPs::test
