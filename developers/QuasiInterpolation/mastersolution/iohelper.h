/**
 * @file iohelper.h
 * @brief NPDE exam TEMPLATE HEADER FILE
 * @author Oliver Rietmann
 * @date 17.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef IOHELPER_H_
#define IOHELPER_H_

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <string>

namespace QuasiInterpolation {

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

void writeCSV(const Eigen::VectorXd& meshwidth, const Eigen::VectorXd& l2_error,
              const Eigen::VectorXd& h1_error, const std::string& filename) {
  std::ofstream file;
  file.open(filename);
  file << meshwidth.transpose().format(CSVFormat) << std::endl;
  file << l2_error.transpose().format(CSVFormat) << std::endl;
  file << h1_error.transpose().format(CSVFormat) << std::endl;
  file.close();
  std::cout << "Generated " + filename << std::endl;
}

void printError(const Eigen::VectorXd& meshwidth,
                const Eigen::VectorXd& l2_error,
                const Eigen::VectorXd& h1_error, const std::string& title) {
  std::cout << title << std::endl;
  std::cout << "meshwidth: " << meshwidth.transpose().format(CSVFormat)
            << std::endl;
  std::cout << "L2-error:  " << l2_error.transpose().format(CSVFormat)
            << std::endl;
  std::cout << "H1-error:  " << h1_error.transpose().format(CSVFormat)
            << std::endl;
}

}  // namespace QuasiInterpolation

#endif  // #ifndef IOHELPER_H_
