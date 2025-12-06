/**
 * @file fvmshallowwater.h
 * @brief NPDE homework FVMShallowWater code
 * @author Rafl Hiptmair
 * @date June 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef FVMShallowWater_H_
#define FVMShallowWater_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace FVMShallowWater {

/** @brief Check for valid 2-shock solution */
bool isPhysicalTwoShockSolution(Eigen::Vector2d ul, Eigen::Vector2d us,
                                Eigen::Vector2d ur);

bool checkSWERHJC(Eigen::Vector2d ul, Eigen::Vector2d ur, double* speed);

bool checkSWEPhysicalShock(Eigen::Vector2d ul, Eigen::Vector2d ur);

/* SAM_LISTING_BEGIN_1 */
inline Eigen::Vector2d sweFluxFunction(Eigen::Vector2d u) {
  assert((u[0] > 0.0) && "Height must be positive!");
  return {u[1], (u[1] * u[1]) / u[0] + 0.5 * u[0] * u[0]};
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
inline std::pair<double, double> sweLambdas(Eigen::Vector2d u) {
  assert((u[0] > 0.0) && "Height must be positive!");
  const double sqrtu1 = std::sqrt(u[0]);
  const double quot = u[1] / u[0];
  return {quot - sqrtu1, quot + sqrtu1};
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
inline Eigen::Matrix2d sweEV(Eigen::Vector2d u) {
  assert((u[0] > 0.0) && "Height must be positive!");
  Eigen::Matrix2d EVs;
  EVs(0, 0) = 1;
  EVs(1, 0) = u[1] / u[0] - std::sqrt(u[0]);
  EVs(0, 1) = 1;
  EVs(1, 1) = u[1] / u[0] + std::sqrt(u[0]);
  return EVs;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
inline Eigen::Matrix2d sweJacobian(Eigen::Vector2d u) {
  assert((u[0] > 0.0) && "Height must be positive!");
  const double quot = u[1] / u[0];
  Eigen::Matrix2d DF;
  DF(0, 0) = 0.0;
  DF(1, 0) = -quot * quot + u[0];
  DF(0, 1) = 1.0;
  DF(1, 1) = 2 * quot;
  return DF;
}
/* SAM_LISTING_END_4 */

/** @brief Coputation of Roe average for SWE and of associated eigenvalues of
 * the Jacobian
 *
 * @param v "left" state
 * @param w "right" state
 * @return triple: (Roe average, lambda1, lambda2)
 * Note that all states are in conservative variables.
 */
std::tuple<Eigen::Vector2d, double, double> RoeAvgSWE(Eigen::Vector2d v,
                                                      Eigen::Vector2d w);

/** @brief HLLE numerical flux for shallow water equations
 *
 * @param v,w the two adjacent cell states
 */
Eigen::Vector2d numfluxHLLESWE(Eigen::Vector2d v, Eigen::Vector2d w);

/** @brief Lax-Friedrichs numerical flux for shallow water equations
 *
 * @param v,w the two adjacent cell states
 */
Eigen::Vector2d numfluxLFSWE(Eigen::Vector2d v, Eigen::Vector2d w);

/** @brief Simple fuly discrete finite volume method for scalar conservation law
 * with explicit Euler timestepping
 *
 * @tparam STATE type for state vectors
 * @tparam TIMESTEPPER functor type for time-local timestep size selection
 * @tparam NUMFLUX functor for numerical flux
 * @tparam RECORDER type for capturing states during timestepping
 *
 * @param a left boundary of computational domain
 * @param b right bounfdary of computational domain
 * @param T final time

 * @param ts_conmtrol  functor for time-local timestep size selection
 * @param u initial data
 * @param NF numerical flux functor
 * @param rec tracker for the FVM fully discrete evolution
 */
/* SAM_LISTING_BEGIN_5 */
template <
    typename STATE, typename TIMESTEPPER, typename NUMFLUX,
    typename RECORDER = std::function<void(double, const std::vector<STATE>&)>>
void FVMEvlGeneric(
    double a, double b, double T, TIMESTEPPER&& ts_control,
    std::vector<STATE>& u, NUMFLUX&& NF,
    RECORDER&& rec = [](double /*t*/, const std::vector<STATE>& /*u*/) -> void {
    }) {
  size_t N = u.size();  // number of cells = number of (spatial) unknowns
  double h = (b - a) / N;
  // Main timestepping loop
  std::vector<STATE> u_tmp{u};
  rec(0.0, u);
  double t = 0.0;
  do {
    u_tmp.swap(u);
    double tau = std::min(ts_control(h, u_tmp), T - t);  // Local timestep size
    // Constant continuation on the left
    STATE F_minus;
    // Implements discrete evolution \lref{eq:Cheeul}
    STATE F_plus = NF(u_tmp[0], u_tmp[0]);
    for (int j = 0; j < N - 1; ++j) {
      F_minus = F_plus;
      F_plus = NF(u_tmp[j], u_tmp[j + 1]);
      u[j] = u_tmp[j] - tau / h * (F_plus - F_minus);
    }
    F_minus = F_plus;
    // Constant continuation on the right
    F_plus = NF(u_tmp[N - 1], u_tmp[N - 1]);
    u[N - 1] = u_tmp[N - 1] - tau / h * (F_plus - F_minus);
    rec(t, u);
    t += tau;
  } while (t < T);
}
/* SAM_LISTING_END_5 */

/** @brief Finite volume solver for 1D shallow water equations on [-1,1]
 *
 * @tparam U0FUNCTOR std::function<Eigen::Vector2d(double)>
 * @tparam NUMFLUX std::function<Eigen::Vector2d(Eigen::Vector2d,
 * Eigen::Vector2d>
 *
 * @double T final time
 * @param N number of mesh cells in [-1,1]
 * @param NF functor object providing numerical flux function
 * @param u0 functor for initial data
 * @param output stream for data
 */

/* SAM_LISTING_BEGIN_6 */
template <typename U0FUNCTOR, typename NUMFLUX>
std::vector<Eigen::Vector2d> solveSWE(double T, unsigned int N, NUMFLUX&& NF,
                                      U0FUNCTOR&& u0,
                                      std::ostream* out = nullptr) {
  // Initialize initial data for N cell values
  std::vector<Eigen::Vector2d> u0_vec{N};
  // Sample initial values on equidistant mesh
  const double h = 2.0 / N;
  double x = -1.0 + h / 2;
  for (int k = 0; k < N; ++k, x += h) {
    u0_vec[k] = u0(x);
  }
  // Functor for timestep control
  auto ts_ctrl = [](double h,
                    const std::vector<Eigen::Vector2d>& u_vec) -> double {
    double s_max = 0.0;
    double l1;
    double l2;
    // Determine maximal characteristic speed
    for (const Eigen::Vector2d& u : u_vec) {
      std::tie(l1, l2) = sweLambdas(u);
      s_max = std::max({s_max, std::abs(l1), std::abs(l2)});
    }
    return h / s_max;
  };
  std::vector<std::vector<Eigen::Vector2d>> data;
  std::vector<double> times;
  FVMEvlGeneric(
      -1.0, 1.0, T, ts_ctrl, u0_vec, NF,
      [&data, &times](double t, const std::vector<Eigen::Vector2d>& u) -> void {
        times.push_back(t);
        data.push_back(u);
      });
  if (out) {
    (*out) << "times = [";
    for (auto t : times) {
      (*out) << t << ' ';
    }
    (*out) << "];\n u1 = [ ...\n";
    for (auto u_vec : data) {
      for (const Eigen::Vector2d& u : u_vec) {
        (*out) << u[0] << ' ';
      }
      (*out) << "; ...\n";
    }
    (*out) << "]; \n u2 = [ ...\n";
    for (auto u_vec : data) {
      for (const Eigen::Vector2d& u : u_vec) {
        (*out) << u[1] << ' ';
      }
      (*out) << "; ...\n";
    }
    (*out) << "];\n";
  }
  return data.back();
}
/* SAM_LISTING_END_6 */

}  // namespace FVMShallowWater

#endif
