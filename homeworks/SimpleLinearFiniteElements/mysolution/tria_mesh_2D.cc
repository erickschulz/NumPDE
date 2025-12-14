#include "tria_mesh_2D.h"

#include <cassert>
#include <fstream>
#include <iostream>

namespace SimpleLinearFiniteElements {

template <typename Derived>
std::istream& operator>>(std::istream& is, Eigen::MatrixBase<Derived>& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      is >> matrix(i, j);
    }
  }
  return is;
}

/**
 * @brief read in mesh file into a costum data structure for this exercise
 * @param filename name of the mesh file
 */
TriaMesh2D::TriaMesh2D(std::string filename) {
  std::ifstream file(filename, std::ifstream::in);

  if (file.is_open()) {
    // read vertices
    int n_vertices;
    std::string vertices_label;
    file >> n_vertices >> vertices_label;
    if (vertices_label != "Vertices") {
      file.close();
      throw std::runtime_error(
          "Keyword 'Vertices' not found. Wrong file format.");
    }
    _nodecoords = Eigen::MatrixXd(n_vertices, 2);
    file >> _nodecoords;

    // read elements
    int n_elements;
    std::string elements_label;
    file >> n_elements >> elements_label;
    if (elements_label != "Elements") {
      file.close();
      throw std::runtime_error(
          "Keyword 'Elements' not found. Wrong file format.");
    }
    _elements = Eigen::MatrixXi(n_elements, 3);
    file >> _elements;

    file.close();
  } else {
    throw std::runtime_error("Error when opening file '" + filename + "'.");
  }
}

TriGeo_t TriaMesh2D::getVtCoords(int i) const {
  TriGeo_t triangle;
  for (int k = 0; k < 3; ++k) {
    triangle.col(k) = _nodecoords.row(_elements(i, k));
  }
  return triangle;
}

/**
 * @brief Saves a 3D mesh file for plotting function on the 2D mesh.
 * The new z component contains the values of the function on the 2D mesh.
 * @param filename Output file name of the new 3D mesh
 * @param z vector of z values, in correct order
 */
void TriaMesh2D::SaveMesh3D(std::string filename,
                            const Eigen::VectorXd& z) const {
  int n_vertices = _nodecoords.rows();
  int n_elements = _elements.rows();

  Eigen::MatrixXd new_vertices(n_vertices, 3);
  new_vertices.leftCols<2>() = _nodecoords;
  new_vertices.col(2) = z;

  std::ofstream file(filename);
  if (file.is_open()) {
    file << n_vertices << " Vertices" << std::endl;
    file << new_vertices << std::endl;

    file << n_elements << " Elements" << std::endl;
    file << _elements;

    file.close();
  } else {
    throw std::runtime_error("Error when opening file '" + filename + "'.");
  }
}

}  // namespace SimpleLinearFiniteElements
