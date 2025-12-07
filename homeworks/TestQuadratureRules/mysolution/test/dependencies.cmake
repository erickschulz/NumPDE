# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/testquadrules_test.cc
)

set(LIBRARIES Eigen3::Eigen GTest::gtest_main LF_ALL)
