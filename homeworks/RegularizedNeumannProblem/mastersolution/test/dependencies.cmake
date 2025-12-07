# Dependencies of mastersolution tests:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/regularizedneumannproblem_test.cc
)

set(LIBRARIES Eigen3::Eigen GTest::gtest_main LF_ALL)
