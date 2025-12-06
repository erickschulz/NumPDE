#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Add your custom dependencies here:
#endif

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/stableevaluationatapoint_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.fe
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.quad
  LF::lf.uscalfe
)

