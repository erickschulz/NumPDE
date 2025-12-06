#ifndef SYMPLECTIC_ASSEMBLE_HPP
#define SYMPLECTIC_ASSEMBLE_HPP
/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

// common includes
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace SymplecticTimesteppingWaves {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNC_ALPHA, typename FUNC_GAMMA, typename FUNC_BETA>
Eigen::SparseMatrix<double> assembleGalerkinMatrix(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNC_ALPHA&& alpha, FUNC_GAMMA&& gamma, FUNC_BETA&& beta) {
  Eigen::SparseMatrix<double> galMat;

  //====================
  // Your code goes here
  //====================
  return galMat;
}
/* SAM_LISTING_END_1 */

}  // namespace SymplecticTimesteppingWaves

#endif
