/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "simplelinearfiniteelements.h"

namespace SimpleLinearFiniteElements {

/* SAM_LISTING_BEGIN_9 */
double getArea(const TriGeo_t& triangle) {
  return std::abs(
      0.5 *
      ((triangle(0, 1) - triangle(0, 0)) * (triangle(1, 2) - triangle(1, 1)) -
       (triangle(0, 2) - triangle(0, 1)) * (triangle(1, 1) - triangle(1, 0))));
}
/* SAM_LISTING_END_9 */

/**
 * @brief Compute the barycentric coordinate function gradients
 * @param vertices The vertices of the triangle
 * @returns A matrix with the gradients in its columns
 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& vertices) {
  Eigen::Matrix3d X;
  // Argument \texttt{vertices} passes the vertex positions of the triangle
  // as the \textbf{columns} of a $2\times 3$-matrix, see
  // \cref{cpp:getVtCoords}. The function returns the components of the
  // gradients as the \textbf{columns} of a $2\times 3$-matrix.

  // Computation based on \eqref{eq:lambdalse}, solving for the
  // coefficints of the barycentric oordinate functions.
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

/**
 *  @brief Computation of element mass matrix on planar triangle
 *  @param V 2x3 matrix of vertex coordinates
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix3d ElementMatrix_Mass_LFE(const TriGeo_t& V) {
  Eigen::Matrix3d element_matrix;
  element_matrix << 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0;
  element_matrix *= getArea(V) / 12.0;
  return element_matrix;
}
/* SAM_LISTING_END_1 */

/**
 *  @brief Computation of Element Matrix for the Laplacian
 *  @param V The vertices of the triangle
 */
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& V) {
  // Argument \texttt{V} same as \texttt{vertices} in \cref{cpp:gradbarycords}.
  // The function returns the $3\times 3$ element matrix as a fixed size
  // \eigen matix.

  Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(V);
  // compute inner products of gradients through matrix multiplication
  return getArea(V) * X.transpose() * X;
}

/**
 *  @brief Computation of full Element Matrix
 *  @param V The vertices of the triangle
 */
Eigen::Matrix3d getElementMatrix(const TriGeo_t& V) {
  return ElementMatrix_Lapl_LFE(V) + ElementMatrix_Mass_LFE(V);
}

/**
 * @brief Compute the local Element Load Vector
 * @param V Vertex coordinates of the triangle
 * @param FHandle Load function f
 * @returns The local Element Load Vector
 */
Eigen::Vector3d localLoadLFE(const TriGeo_t& V, const FHandle_t& FHandle) {
  Eigen::Vector3d philoc = Eigen::Vector3d::Zero();
  // Evaluate source function for ertex locations
  for (int i = 0; i < 3; ++i) {
    philoc(i) = FHandle(V.col(i));
  }
  philoc *= getArea(V) / 3;
  return philoc;
}

/**
 * @brief Assembles the Galerkin Matrix
 * @param mesh the mesh to use
 * @param getElementMatrix Element Matrix
 * @return Galerkin Matrix
 */
Eigen::SparseMatrix<double> assembleGalMatLFE(
    const TriaMesh2D& Mesh, const LocalMatrixHandle_t& getElementMatrix) {
  // obtain the number of vertices
  int N = Mesh._nodecoords.rows();
  // obtain the number of elements/cells
  int M = Mesh._elements.rows();
  std::vector<Eigen::Triplet<double>> triplets;
  // loop over elements and add local contributions
  for (int i = 0; i < M; i++) {
    // get local$\to$global index mapping for current element, \emph{cf.}
    // \lref{eq:idxdef}
    const TriGeo_t triangle = Mesh.getVtCoords(i);
    // compute element contributions
    Eigen::Matrix3d Ak = getElementMatrix(triangle);
    // build triplets from contributions
    Eigen::Vector3i element = Mesh._elements.row(i);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        triplets.emplace_back(element(j), element(k), Ak(j, k));
      }
    }
  }
  // build sparse matrix from triplets
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}

/**
 * @brief assemLoad_LFE Assembles the Load Vector
 * @param mesh the mesh to use
 * @param getElementVector Element Vector factory
 * @param FHandle function handle for f
 * @return assembled load vector
 */
Eigen::VectorXd assemLoad_LFE(const TriaMesh2D& Mesh,
                              const LocalVectorHandle_t& getElementVector,
                              const FHandle_t& FHandle) {
  // obtain the number of triangles
  int M = Mesh._elements.rows();

  // obtain the number of vertices
  int N = Mesh._nodecoords.rows();
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N);

  // loop over all triangles
  for (int i = 0; i < M; i++) {
    const TriGeo_t Vertices = Mesh.getVtCoords(i);
    // Compute the element right hand side vector
    const Eigen::Vector3d philoc = getElementVector(Vertices, FHandle);
    // Add contibutions to global load vector
    const Eigen::Vector3i dofhk = Mesh._elements.row(i);
    for (int j = 0; j < 3; ++j) {
      phi(dofhk(j)) += philoc(j);
    }
  }

  return phi;
}

/**
 * @brief H1Serror Computes the H^1 error between the approximate solution and
 *                the exact solution
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact the exact gradient of the solution
 * @return the H^1 difference
 */
/* SAM_LISTING_BEGIN_3 */
double H1Serror(
    const TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d&)> exact) {
  double H1Serror_squared = 0.0;

  // loop over all triangles
  for (int i = 0; i < mesh._elements.rows(); ++i) {
    const TriGeo_t triangle = mesh.getVtCoords(i);

    // loop over all three vertices of the triangle
    Eigen::Vector3d values_at_vertices;
    for (int k = 0; k < 3; ++k) {
      values_at_vertices(k) = uFEM[mesh._elements(i, k)];
    }

    // gradient of FEM approximation (same for all 3 vertices!)
    Eigen::Vector2d gradient_FEM =
        gradbarycoordinates(triangle) * values_at_vertices;

    // loop over all three vertices of the triangle
    Eigen::Vector3d error_at_vertices;
    for (int k = 0; k < 3; ++k) {
      Eigen::Vector2d gradient_exact = exact(triangle.col(k));
      error_at_vertices(k) = (gradient_FEM - gradient_exact).squaredNorm();
    }

    // Add squared error per triangle
    H1Serror_squared += getArea(triangle) / 3.0 * error_at_vertices.sum();
  }

  return std::sqrt(H1Serror_squared);
}
/* SAM_LISTING_END_3 */

/**
 * @brief solves system and prints H1-semierror, L2 error, the mesh and a
 * surface plot
 * @param mesh discretisation of the computational domain
 * @returns A tuple of the solution, the L2-Error and the H1-Error
 */
/* SAM_LISTING_BEGIN_4 */
std::tuple<Eigen::VectorXd, double, double> Solve(
    const SimpleLinearFiniteElements::TriaMesh2D& mesh) {
  const double pi = 3.1415926535897;

  // define the source function f
  auto f = [pi](const Eigen::Vector2d& x) {
    return (1.0 + 8.0 * pi * pi) * std::cos(2.0 * pi * x(0)) *
           std::cos(2.0 * pi * x(1));
  };
  // the exact solution of the linear variational problem
  auto uExact = [pi](const Eigen::Vector2d& x) {
    return std::cos(2 * pi * x(0)) * std::cos(2 * pi * x(1));
  };

  Eigen::VectorXd U;
  double l2error;
  double h1error;

  // the gradient of uExact that can be easily analytically computed
  auto gradUExact = [pi](const Eigen::Vector2d& x) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x(0)) * std::cos(2 * pi * x(1)),
        -2 * pi * std::cos(2 * pi * x(0)) * std::sin(2 * pi * x(1));
    return gradient;
  };
  // Compute the galerkin matrix A and load vector L
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::assembleGalMatLFE(
      mesh, SimpleLinearFiniteElements::getElementMatrix);

  Eigen::VectorXd L =
      SimpleLinearFiniteElements::assemLoad_LFE(mesh, localLoadLFE, f);

  // solve the LSE using the sparse LU solver of Eigen
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  U = solver.solve(L);

  // compute both the L2 and the H1s error
  l2error = SimpleLinearFiniteElements::L2Error(mesh, U, uExact);
  h1error = SimpleLinearFiniteElements::H1Serror(mesh, U, gradUExact);
  return std::make_tuple(U, l2error, h1error);
}
/* SAM_LISTING_END_4 */

}  // namespace SimpleLinearFiniteElements
