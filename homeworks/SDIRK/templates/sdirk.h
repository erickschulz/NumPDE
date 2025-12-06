#ifndef SDIRK_H_
#define SDIRK_H_

/**
 * @file sdirk.h
 * @brief NPDE homework SDIRK code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <vector>

namespace SDIRK {

// Compute one step of  the IVP y'' + y' + y = 0 using a SDIRK method
Eigen::Vector2d SdirkStep(const Eigen::Vector2d& z0, double h, double gamma);

// Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK
std::vector<Eigen::Vector2d> SdirkSolve(const Eigen::Vector2d& z0,
                                        unsigned int M, double T, double gamma);

double CvgSDIRK();

}  // namespace SDIRK

#endif  // #define SDIRK_H_
