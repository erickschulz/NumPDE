/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */
#ifndef SLFEM2D
#define SLFEM2D

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <exception>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

/**
 * @brief The TriaMesh2D struct describes the triangulation in the form of a
 * _nodecoords matrix that holds the coordinates of vertex i in its i-th row and
 * an _elements matrix that holds the three indices of the vertices forming a
 * triangle in each row
 */
/* SAM_LISTING_BEGIN_1 */
// Matrix containing vertex coordinates of a triangle
using TriGeo_t = Eigen::Matrix<double, 2, 3>;
struct TriaMesh2D {
  // Constructor: reads mesh data from file, whose name is passed
  TriaMesh2D(std::string filename);  // \Label[line]{tm:cs}
  virtual ~TriaMesh2D(void) {}
  // Retrieve coordinates of vertices of a triangles as columns
  // of a fixed-size 2x3 matrix
  TriGeo_t getVtCoords(std::size_t) const;
  // Data members describing geometry and topolgy
  Eigen::Matrix<double, Eigen::Dynamic, 2> _nodecoords;
  Eigen::Matrix<int, Eigen::Dynamic, 3> _elements;
};
/* SAM_LISTING_END_1 */

// Signature of a function computing the element matrix for a triangular cell
// and piecewise linear Lagrangian finite elements
typedef std::function<Eigen::Matrix3d(const TriGeo_t&)> LocalMatrixHandle_t;
// Data type for a COO matrix: sequence of triplets
typedef std::vector<Eigen::Triplet<double>> Triplet_t;

/**
 * @brief The MatrixAssembler class assembles a Galerkin matrix given a
 * function for computing local matrices from TriGeo_t triangles.
 */
class MatrixAssembler {
 public:
  // Constructor: stores element matrix assembler
  MatrixAssembler(LocalMatrixHandle_t getElementMatrix)
      : localMatrixHandle(std::move(getElementMatrix)) {}

  // Assemble the Galerkin matrix for the provided mesh
  Eigen::SparseMatrix<double> Assemble(const TriaMesh2D& mesh);

 private:
  LocalMatrixHandle_t localMatrixHandle;
};

/**
 * @brief The SlowMatrixAssembler class assembles a Galerkin matrix given a
 * function for computing local matrices from TriGeo_t triangles using a less
 * efficient way to build the sparse matrix.
 */
class SlowMatrixAssembler {
 public:
  // Constructor: stores element matrix assembler
  SlowMatrixAssembler(LocalMatrixHandle_t getElementMatrix)
      : localMatrixHandle(std::move(getElementMatrix)) {}

  // Assemble the Galerkin matrix for the provided mesh
  Eigen::SparseMatrix<double> Assemble(const TriaMesh2D& mesh);

 private:
  LocalMatrixHandle_t localMatrixHandle;
};

// Type for a real-valued function on the computational domain
typedef std::function<double(const Eigen::Vector2d&)> FHandle_t;
// Signature of a function computing element vectors
typedef std::function<Eigen::Vector3d(const TriGeo_t&, FHandle_t)>
    LocalVectorHandle_t;

/**
 * @brief The VectorAssembler class assembles a load vector given a
 * function for computing local vectors from TriGeo_t triangles.
 */
class VectorAssembler {
 public:
  // Constructor: stores element matrix assembler and load function
  VectorAssembler(const LocalVectorHandle_t getElementVector,
                  FHandle_t sourceFunction)
      : localVectorHandle(std::move(getElementVector)),
        sourceFunction(std::move(sourceFunction)) {}

  // Assemble the load vector for the provided mesh
  Eigen::VectorXd Assemble(const TriaMesh2D& mesh);

 private:
  LocalVectorHandle_t localVectorHandle;
  FHandle_t sourceFunction;
};

/**
 * @brief The finite element solver solves the problem
 * \int_{\Omega} grad(u)grad(v) + uv dx = \int_{\Omega} fv dx
 * using MatrixAssembler and VectorAssembler. If the inefficient_flag is
 * set SlowMatrixAssembler is used to demonstrate the negative effect of not
 * building the Galerkin matrix from triplets.
 */
class FESolver {
 public:
  // Constructor: stores source function f and a flag indicating
  // if the slow Galerkin matrix assembler should be used
  FESolver(const FHandle_t& sourceFunction, int inefficient_flag = 0)
      : sourceFunction(sourceFunction) {
    inefficient = inefficient_flag;
  };

  // Solve the discretized system
  Eigen::VectorXd Solve(const TriaMesh2D& mesh);

 private:
  const FHandle_t& sourceFunction;
  int inefficient;
};
#endif
