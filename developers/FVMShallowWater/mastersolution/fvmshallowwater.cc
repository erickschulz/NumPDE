/**
 * @file fvmshallowwater.cc
 * @brief NPDE homework FVMShallowWater code
 * @author Ralf Hiptmair
 * @date June 2024
 * @copyright Developed at SAM, ETH Zurich
 */

#include "fvmshallowwater.h"

#include <cmath>
#include <limits>

namespace FVMShallowWater {

#if SOLUTION
/* SAM_LISTING_BEGIN_2 */
std::tuple<Eigen::Vector2d, double, double> RoeAvgSWE(Eigen::Vector2d v,
                                                      Eigen::Vector2d w) {
  assert((v[0] > 0.0) && "v-height must be positive!");
  assert((w[0] > 0.0) && "w-height must be positive!");
  // Roe average state for SWE
  Eigen::Vector2d u_bar;
  // Implements \prbeqref{eq:ras} to compute $\overline{\Vu}$
  u_bar[0] = 0.5 * (v[0] + w[0]);  // Arithmetic average of first component
  const double sv0 = std::sqrt(v[0]);
  const double sw0 = std::sqrt(w[0]);
  // Row average in physical velocity
  const double v_avg = (sv0 * v[1] / v[0] + sw0 * w[1] / w[0]) / (sv0 + sw0);
  u_bar[1] = v_avg * u_bar[0];
  const std::pair<double, double> lambdas{sweLambdas(u_bar)};
  return {u_bar, lambdas.first, lambdas.second};
}
/* SAM_LISTING_END_2 */
#endif

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d numfluxLFSWE(Eigen::Vector2d v, Eigen::Vector2d w) {
  assert((v[0] > 0.0) && "v-height must be positive!");
  assert((w[0] > 0.0) && "w-height must be positive!");
  Eigen::Vector2d nfLF{};
#if SOLUTION
  const auto& F = &sweFluxFunction;
  // Compute Row average and asscoiated characteristic speeds
  const auto [u_bar, l1, l2] = RoeAvgSWE(v, w);
  // Compute "modulus" of Roe matrix
  const Eigen::Matrix2d EV_bar = sweEV(u_bar);
  const Eigen::Matrix2d D_abs =
      (Eigen::Matrix2d() << std::abs(l1), 0.0, 0.0, std::abs(l2)).finished();
  const Eigen::Matrix2d AR_abs = (EV_bar * D_abs) * EV_bar.inverse();
  // Implment \lref{eq:lfgsys}
  nfLF = 0.5 * (F(v) + F(w)) - 0.5 * AR_abs * (w - v);
#else
  /* **************************************************
     Your code here
     ************************************************** */
#endif
  return nfLF;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
Eigen::Vector2d numfluxHLLESWE(Eigen::Vector2d v, Eigen::Vector2d w) {
  assert((v[0] > 0.0) && "v-height must be positive!");
  assert((w[0] > 0.0) && "w-height must be positive!");
  Eigen::Vector2d nfHLLE{};
#if SOLUTION
  const auto& F = &sweFluxFunction;
  // Compute characteristic speeds for edge states
  const auto [l1v, l2v] = sweLambdas(v);
  const auto [l1w, l2w] = sweLambdas(w);
  // Compute Row average and associated characteristic speeds
  auto [u_bar, l1, l2] = RoeAvgSWE(v, w);
  // HLLE edge shock speeds according ro \lref{eq:HLLes}
  const double s_minus = std::min(l1, l1v);
  const double s_plus = std::max(l2, l2w);
  // Compute middle state $\Vu^*$ according to \lref{eq:usteq}
  const Eigen::Vector2d Fv{F(v)};
  const Eigen::Vector2d Fw{F(w)};
  if (s_minus > 0) {
    nfHLLE = Fv;
  } else if (s_plus < 0) {
    nfHLLE = Fw;
  } else {
    const Eigen::Vector2d us =
        (Fw - Fv - s_plus * w + s_minus * v) / (s_minus - s_plus);
    nfHLLE = F(us);
  }
#else
  /* **************************************************
     Your code here
     ************************************************** */
#endif
  return nfHLLE;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
#if SOLUTION
bool checkSWERHJC(Eigen::Vector2d ul, Eigen::Vector2d ur, double* speed) {
  Eigen::Matrix2d M;
  double g = 1;
  M << ur(0) - ul(0), ul(1) - ur(1), ur(1) - ul(1),
      ul(1) * ul(1) / ul(0) + .5 * g * ul(0) * ul(0) -
          (ur(1) * ur(1) / ur(0) + .5 * g * ur(0) * ur(0));
  auto rank = M.jacobiSvd()
                  .setThreshold(10. * std::numeric_limits<double>::epsilon())
                  .rank();
  if (rank < 2) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
    Eigen::MatrixXd M_null_space = lu.kernel();
    *speed = M_null_space(0) / M_null_space(1);
    return true;
  }
  return false;
}
#endif
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
#if SOLUTION
bool checkSWEPhysicalShock(Eigen::Vector2d ul, Eigen::Vector2d ur) {
  double s;
  auto lambda_1 = [](Eigen::Vector2d u) -> double {
    return u(1) / u(0) - std::sqrt(u(0));
  };
  auto lambda_2 = [](Eigen::Vector2d u) -> double {
    return u(1) / u(0) + std::sqrt(u(0));
  };
  if (checkSWERHJC(ul, ur, &s)) {
    // find k
    if (lambda_1(ul) > s && s > lambda_1(ur)) {
      if (lambda_2(ul) > s && lambda_2(ur) > s) return true;
    } else if (lambda_2(ul) > s && s > lambda_2(ur)) {
      if (lambda_1(ul) < s && lambda_1(ur) < s) return true;
    }
  }
  return false;
}
#endif
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
bool isPhysicalTwoShockSolution(Eigen::Vector2d ul, Eigen::Vector2d us,
                                Eigen::Vector2d ur) {
#if SOLUTION
  // Check if the states are admissible
  if (ul(0) <= 0 || us(0) <= 0 || ur(0) <= 0) return false;

  // Check whether the states are equal
  bool ul_eq_us =
      (ul - us).norm() < 10. * std::numeric_limits<double>::epsilon();
  bool us_eq_ur =
      (us - ur).norm() < 10. * std::numeric_limits<double>::epsilon();

  // We check the Rankine-Hugoniot Jump Conditions
  double s1, s2;
  if (!ul_eq_us && !checkSWERHJC(ul, us, &s1)) return false;
  if (!us_eq_ur && !checkSWERHJC(us, ur, &s2)) return false;
  if (!ul_eq_us && !us_eq_ur && s1 > s2) return false;

  // We check the Lax Entropy Condition
  if (!checkSWEPhysicalShock(ul, us)) return false;
  if (!checkSWEPhysicalShock(us, ur)) return false;

    // If all tests passed, we have a physical shock
#endif
  return true;
}
/* SAM_LISTING_END_6 */

}  // namespace FVMShallowWater
