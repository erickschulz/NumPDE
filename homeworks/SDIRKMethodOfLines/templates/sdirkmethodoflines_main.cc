/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
#include <string>

#include "sdirkmethodoflines.h"
#include "sdirkmethodoflines_ode.h"
#include "systemcall.h"

using namespace SDIRKMethodOfLines;

int main(int /*argc*/, char** /*argv*/) {
  /* SDIRK-2 ODE convergence */
  sdirk2ScalarODECvTest();

  std::cout << "\n************* SDIRKMethodOfLines *************" << std::endl;
  /* Solving the parabolic temperature convection cooling problem */
  // Create a Lehrfem++ square tensor product mesh
  // Obtain mesh factory
  /* std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
       std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
   // Triangular tensor product mesh
   lf::mesh::hybrid2d::TPTriagMeshBuilder builder(mesh_factory_ptr);
   // Set mesh parameters following the Builder pattern
   // Domain is the unit square
   builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
       .setTopRightCorner(Eigen::Vector2d{1, 1})
       .setNoXCells(100)
       .setNoYCells(100);
   std::shared_ptr<lf::mesh::Mesh> mesh_p{builder.Build()}; */

  //====================
  // Your code goes here
  //====================

//====================
// Your code goes here
//====================
  return 0;
}
