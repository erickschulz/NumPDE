/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves.h"

#include <fstream>

namespace SymplecticTimesteppingWaves {

/* SAM_LISTING_BEGIN_7 */
void wavePropSimulation(unsigned int m) {
#if SOLUTION
  double T = 10.0;  // final time

  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), "meshes/hex4.msh");
  auto mesh_p = reader.mesh();

  // Finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Building initial condition vector
  Eigen::VectorXd u0_vec(N_dofs);
  Eigen::VectorXd v0_vec(N_dofs);
  u0_vec.setOnes();
  v0_vec.setZero();

  // Creating coefficient function
  auto c = [](Eigen::Vector2d x) -> double { return 1.0 + x.dot(x); };

  // Solve the initial value problem
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solution_pair =
      solvewave(fe_space, c, u0_vec, v0_vec, T, m);
  Eigen::VectorXd discrete_wave_sol = solution_pair.first;
  LF_ASSERT_MSG(
      discrete_wave_sol.size() == N_dofs,
      "Size of discrete solution and dimension of FE space mismatch.");
  Eigen::VectorXd energies = solution_pair.second;
  LF_ASSERT_MSG(energies.size() == m + 1, "Wrong number of energie values.");

  // Define output file format for the energies
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::string errors_file_name = "waves_energies.csv";
  std::ofstream file(errors_file_name.c_str());
  if (file.is_open()) {
    file << energies.format(CSVFormat);
  }

  std::cout << ">>The energies were written to:" << std::endl;
  std::cout << "\t waves_energies.csv" << std::endl;
  /* SAM_LISTING_BEGIN_2 */
  // Output results for the wave function to vtk file
  lf::io::VtkWriter vtk_writer(mesh_p, "discrete_wave_sol.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        discrete_wave_sol[global_idx];
  };
  vtk_writer.WritePointData("discrete_wave_sol", *nodal_data);
  /* SAM_LISTING_END_2 */
  std::cout << ">>The output discrete_wave_sol was written to:" << std::endl;
  std::cout << "\t discrete_wave_sol.vtk" << std::endl;
#else
  //====================
  // Your code goes here
  //====================
#endif
}

/* SAM_LISTING_END_7 */

void progress_bar::write(double fraction) {
  // clamp fraction to valid range [0,1]
  if (fraction < 0)
    fraction = 0;
  else if (fraction > 1)
    fraction = 1;

  auto width = bar_width - message.size();
  auto offset = bar_width - static_cast<unsigned>(width * fraction);

  os << '\r' << message;
  os.write(full_bar.data() + offset, width);
  os << " [" << std::setw(3) << static_cast<int>(100 * fraction) << "%] "
     << std::flush;
}

/* SAM_LISTING_BEGIN_5 */
double testStab() {
  double maxUniformTimestep;
#if SOLUTION
  std::cout << "\n*********************** STABILITY EXPERIMENT "
               "***********************"
            << std::endl;
  unsigned int m_upper = 1100;
  unsigned int m_lower = 1000;
  unsigned int m = m_lower + (m_upper - m_lower) / 2;
  while (m_upper - m_lower > 1) {
    std::cout << "\nTesting symplectic method with m = " << m << std::endl;
    // Catch exception thrown in case of blow-up
    try {
      wavePropSimulation(m);
    } catch (const char* msg) {
      // Blow-up detected!
      std::cout << "Energy blows up!" << std::endl;
      m_lower = m;
      m = m_lower + (m_upper - m_lower) / 2;
      continue;
    }
    // No blow-up
    std::cout << "Energy is conserved" << std::endl;
    m_upper = m;
    m = m_lower + (m_upper - m_lower) / 2;
  }
  maxUniformTimestep = 10.0 / m_upper;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return maxUniformTimestep;
}
/* SAM_LISTING_END_5 */

}  // namespace SymplecticTimesteppingWaves
