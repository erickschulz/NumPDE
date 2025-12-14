#ifndef NLMATODE_H_
#define NLMATODE_H_

#include <Eigen/Core>

namespace NLMatODE {

//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
Eigen::MatrixXd matode(const Eigen::MatrixXd& Y0, double T);

//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
bool checkinvariant(const Eigen::MatrixXd& M, double T);

double cvgDiscreteGradientMethod();

}  // namespace NLMatODE

#endif  // #define NLMATODE_H_
