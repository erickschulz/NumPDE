/**
 * @ file transformedconslaw.h
 * @ brief NPDE homework about conservation law with non-linear density
 * @ author Ralf Hiptmair, Oliver Rietmann
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <utility>

namespace TRFCL {

/**
 * @brief Class describing a Cauchy problem with non-linear density
 */
class NonStdCauchyProblemCL {
 public:
  NonStdCauchyProblemCL() = default;
  virtual ~NonStdCauchyProblemCL() = default;
  // Function rho
  double rho(double z) const;
  // Derivative of function rho
  double drho(double z) const;
  // Function g
  double g(double z) const;
  // Derivative of function g
  double dg(double z) const;
  // Finite interval containing "interesting" parts of solution
  std::pair<double, double> domain() const;
  // Final time
  double T() const;
  // Initial data
  double z0(double x) const;
};

/**
 * @brief Newton's method for the inversion of rho
 *
 * Approximately solves rho(z) = u employing Newton's method with initial guess.
 *
 * @tparam RHOFUNCTOR provides rho(z), std::function<double(double)>
 * @tparam DRHOFUNCTOR provides rho'(z), std::function<double(double)>
 * @param u value for which the the inverse of rho should be computed
 * @param z0 initial guess for Newton's method
 * @param rho functor for z -> rho(z)
 * @param rhod functor providing the derivative of rho
 * @param atol absolute tolerance
 * @param rtol relative tolerance

 WRONG IMPLEMENTATION OF TEARMINATION CRITERION.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename RHOFUNCTOR, typename DRHOFUNCTOR>
double rhoInverse(double u, double z0, RHOFUNCTOR&& rho, DRHOFUNCTOR&& drho,
                  double atol = 1.0E-10, double rtol = 1.0E-5) {
#if SOLUTION
  // Reduce to problem of finding zeros: phi(z)=0 <==> rho(z)=u
  auto F = [rho, u](double z) -> double { return rho(z) - u; };
  auto dF = [drho](double z) -> double { return drho(z); };

  // Newton iteration for finding a zero of F
  double s;
  do {
    s = F(z0) / dF(z0);
    z0 = z0 - s;
  } while (std::abs(s) > rtol * std::abs(z0) && std::abs(s) > atol);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return z0;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <class CAUCHYPROBLEM>
Eigen::VectorXd semiDiscreteRhs(const Eigen::VectorXd& mu,
                                const Eigen::VectorXd& zeta,
                                CAUCHYPROBLEM prb) {
  int N = mu.size();
  Eigen::VectorXd rhs(N);

#if SOLUTION
  // Define the (simplified) Rusanov 2-point numerical flux
  auto F = [&mu, &zeta, &prb](int i, int j) {
    double v = mu[i];
    double w = mu[j];
    double fv = prb.g(zeta[i]);
    double fw = prb.g(zeta[j]);
    double dfv = prb.dg(zeta[i]) / prb.drho(zeta[i]);
    double dfw = prb.dg(zeta[j]) / prb.drho(zeta[j]);
    double max_abs_df = std::max(std::abs(dfv), std::abs(dfw));
    return 0.5 * (fv + fw) - 0.5 * (w - v) * max_abs_df;
  };

  // Compute (-h) * [right-hand side]
  rhs[0] = F(0, 1) - F(0, 0);
  for (int j = 1; j < N - 1; ++j) {
    rhs[j] = F(j, j + 1) - F(j - 1, j);
  }
  rhs[N - 1] = F(N - 1, N - 1) - F(N - 2, N - 1);

  // Compute cell spacing h
  std::pair<double, double> limits = prb.domain();
  double h = (limits.second - limits.first) / (N - 1);

  // Include factor -1.0 / h
  rhs = (-1.0 / h) * rhs;
#else
  //====================
  // Your code goes here
  //====================
#endif

  return rhs;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <class CAUCHYPROBLEM, typename RECORDER = std::function<
                                   void(double, const Eigen::VectorXd&)>>
Eigen::VectorXd solveCauchyPrb(
    unsigned int M, unsigned int N, CAUCHYPROBLEM prb,
    RECORDER&& rec = [](double /*time*/, const Eigen::VectorXd&
                        /*zstate*/) -> void {}) {
  // Get inital data for zeta
  std::pair<double, double> limits = prb.domain();
  Eigen::VectorXd x =
      Eigen::VectorXd::LinSpaced(N, limits.first, limits.second);
  auto z0 = [&prb](double y) { return prb.z0(y); };
  Eigen::VectorXd zeta = x.unaryExpr(z0);

  // Compute time step
  double T = prb.T();
  double dt = T / M;

  // Wrap the involved functions for simpler use below
  auto rho = [&prb](double v) { return prb.rho(v); };
  auto drho = [&prb](double v) { return prb.drho(v); };
  auto r = [rho, drho](double u, double z) {
    return rhoInverse(u, z, rho, drho);
  };

  // Get inital data for mu and get the time grid
  Eigen::VectorXd mu = zeta.unaryExpr(rho);
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(M + 1, 0.0, T);

  // record solution at initial time
  rec(t[0], zeta);

  for (int k = 0; k < M; ++k) {
#if SOLUTION
    // Compute increment k1
    Eigen::VectorXd k1 = semiDiscreteRhs(mu, zeta, prb);

    // Compute incerment k2
    Eigen::VectorXd mu_temp = mu + dt / 3.0 * k1;
    Eigen::VectorXd zeta_temp = mu_temp.binaryExpr(zeta, r);
    Eigen::VectorXd k2 = semiDiscreteRhs(mu_temp, zeta_temp, prb);

    // Compute increment k3
    mu_temp = mu + dt * 2.0 / 3.0 * k2;
    zeta_temp = mu_temp.binaryExpr(zeta, r);
    Eigen::VectorXd k3 = semiDiscreteRhs(mu_temp, zeta_temp, prb);

    // Update mu and zeta according to the Butcher tableau
    mu = mu + dt * (0.25 * k1 + 0.0 * k2 + 0.75 * k3);
    zeta = mu.binaryExpr(zeta, r);
#else
    //====================
    // Your code goes here
    //====================
#endif

    // record solution after current time step
    rec(t[k + 1], zeta);
  }

  return zeta;
}
/* SAM_LISTING_END_3 */

}  // namespace TRFCL
