/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include <Eigen/Core>
#include <Eigen/LU>
#include <fstream>

/**
 * @brief   Functions for local computations on triangles.
 */

typedef std::function<double(const Eigen::Vector2d&)> FHandle_t;
typedef Eigen::Matrix<double, 2, 3> TriGeo_t;
typedef std::function<Eigen::Matrix3d(const TriGeo_t&)> LocalMatrixHandle_t;

// Compute element stiffness matrix on triangle
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t&);
// Compute element mass matrix on triangle
Eigen::Matrix3d ElementMatrix_Mass_LFE(const TriGeo_t&);
// Compute combined stiffness and mass matrix
Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const TriGeo_t& Vertices);
// Compute load vector on triangle
Eigen::Vector3d localLoadLFE(const TriGeo_t& vertices,
                             const FHandle_t& FHandle);
