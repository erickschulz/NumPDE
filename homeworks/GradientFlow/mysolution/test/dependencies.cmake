# Add your custom dependencies here:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/gradientflow_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
