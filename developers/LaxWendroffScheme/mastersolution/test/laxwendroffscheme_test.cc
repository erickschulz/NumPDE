/**
 * @file laxwendroffscheme_test_mastersolution.cc
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "../laxwendroffscheme.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace LaxWendroffScheme::test {

TEST(LaxWendroffScheme, solveLaxWendroff) {
  double T = 1.0;
  unsigned int M = 20;
  Eigen::VectorXd x = getXValues(T, M);
  auto u_initial = [](double x) { return x >= 0 ? 1.0 : 0.0; };
  Eigen::VectorXd u0 = x.unaryExpr(u_initial);

  // to test
  Eigen::VectorXd sol = solveLaxWendroff(u0, T, M);

  // reference
  Eigen::VectorXd sol_ref(54);
  sol_ref << 0, 0, 0, 4.358531747509734e-20, -2.21981371529156e-17,
      2.272959845702702e-15, -9.766575241780475e-14, 2.961566748871577e-12,
      -6.701834376706102e-11, 1.17155735646535e-09, -1.615773168269057e-08,
      1.77909014168506e-07, -1.571225233065245e-06, 1.110583565942615e-05,
      -6.218648147595973e-05, 0.000269773378056873, -0.0008672298148430562,
      0.001858236660158711, -0.001783081817579445, -0.002663707250640729,
      0.009569537728594221, -0.002562507685016063, -0.02311757114950184,
      0.009271879687815858, 0.04857871487087509, 0.002993598278723432,
      -0.08333469151482514, -0.1046353375394496, -0.0539281836955189,
      0.03081060739055947, 0.1249255418885772, 0.2181465307350215,
      0.3068838998432271, 0.3901614970388319, 0.4679787990791593,
      0.5406783658605172, 0.608699582846934, 0.6724826670156486,
      0.7324329060356791, 0.788910155312715, 0.842232799129924,
      0.8927139771102892, 0.9409345530172758, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  double error = (sol - sol_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(LaxWendroffScheme, numexpLaxWendroffRP) {
  Eigen::VectorXi M(6);
  M << 20, 40, 80, 160, 320, 640;
  Eigen::VectorXd diffL1 = numexpLaxWendroffRP(M);

  Eigen::VectorXd diffL1_ref(6);
  diffL1_ref << 0.1436558625060627, 0.07609298911593251, 0.03988655606579881,
      0.0205247202997978, 0.01047772054580396, 0.005316757349856763;

  double error = (diffL1 - diffL1_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(LaxWendroffScheme, referenceSolution) {
  const double T = 1.0;
  unsigned int M = 20;
  Eigen::VectorXd x = getXValues(T, M);

  // to test
  Eigen::VectorXd u = referenceSolution(x);

  // reference
  Eigen::VectorXd u_ref(54);
  // u_ref << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  // 0, 0, 0, 0, 0, 0, 0, 0, 0.000000035677855, 0.013390319307041,
  // 0.061606854493118, 0.123064156650103, 0.189343264885395, 0.256737080477183,
  // 0.323386304178452, 0.388319216219640, 0.451033788354385, 0.511288550309748,
  // 0.568988602367676, 0.624119567758121, 0.676707015335232, 0.726789360009730,
  // 0.774396804183109, 0.819530198137162, 0.862131833974998, 0.902031702002586,
  // 0.938821690859780, 0.971467676054509, 0.996323269358076, 1.0, 1.0, 1.0;  //
  // old matlab solution
  u_ref << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      8.38050202998275e-294, 4.45548455770779e-197, 6.57428646282557e-119,
      1.78277117388876e-61, 3.40490659115364e-27, 3.50975951036953e-17,
      4.33152157119415e-17, 1.85604084958367e-17, 4.28291449600763e-17,
      8.11754214596344e-18, 1.02530856551663e-17, -6.76779690716347e-11,
      -1.69441336451043e-06, 0.0133914423820019, 0.0616073331998487,
      0.123064472598099, 0.189343522542103, 0.256737313738816,
      0.323386525604616, 0.388319429989982, 0.451033994692818,
      0.511288747244176, 0.568988786481324, 0.62411973449188, 0.676707159047955,
      0.726789473865465, 0.774396879860547, 0.819530225287819,
      0.862131799243635, 0.902031587172126, 0.938821469066465,
      0.971467304312063, 0.996322740594709, 1, 1, 1;

  double error = (u - u_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-5;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(LaxWendroffScheme, numexpLaxWendroffSmoothU0) {
  Eigen::VectorXi M(6);
  M << 20, 40, 80, 160, 320, 640;
  Eigen::VectorXd error = numexpLaxWendroffSmoothU0(M);

  Eigen::VectorXd error_ref(6);
  error_ref << 0.0237401770928193, 0.00666074781328248, 0.00186341093077163,
      0.000503507323563576, 0.000134754219994167, 3.47467216817475e-05;

  double max_diff = (error - error_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-5;
  EXPECT_NEAR(max_diff, 0.0, tol);
}

TEST(LaxWendroffScheme, solveGodunov) {
  double T = 1.0;
  unsigned int M = 20;
  Eigen::VectorXd x = getXValues(T, M);

  constexpr double PI = 3.14159265358979323846;
  auto Square = [](double x) { return x * x; };
  auto smoothU0 = [&Square](double x) {
    return x < 0.0 ? 0.0 : 1.0 < x ? 1.0 : Square(std::sin(0.5 * PI * x));
  };
  Eigen::VectorXd u0 = x.unaryExpr(smoothU0);

  // to test
  Eigen::VectorXd u = solveGodunov(u0, T, M);

  // reference
  Eigen::VectorXd u_ref(54);
  u_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00000449518262136707,
      0.0000667840957603485, 0.000477402028680055, 0.00218883970351049,
      0.00725614988345724, 0.0186633263538622, 0.0392093707991762,
      0.0701841501914256, 0.110918874119619, 0.159370222982678,
      0.213048782942228, 0.269689172597778, 0.327528636281691,
      0.385323806127252, 0.442257037315101, 0.497822413221615,
      0.551728034623663, 0.603822898577608, 0.654045678367657,
      0.702390006164755, 0.748881220936979, 0.793560658950387,
      0.836474589001637, 0.877665727963671, 0.917167938264825,
      0.955044660415018, 0.994503134771148, 1.0, 1.0, 1.0;

  double error = (u - u_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(LaxWendroffScheme, numexpGodunovSmoothU0) {
  Eigen::VectorXi M(6);
  M << 20, 40, 80, 160, 320, 640;
  Eigen::VectorXd error = numexpGodunovSmoothU0(M);

  Eigen::VectorXd error_ref(6);
  error_ref << 0.0682749240658965, 0.0355990327324049, 0.0181554698431239,
      0.00917216574843169, 0.00460908965912379, 0.00231033220753803;

  double max_diff = (error - error_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-5;
  EXPECT_NEAR(max_diff, 0.0, tol);
}

}  // namespace LaxWendroffScheme::test
