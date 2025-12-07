# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/regularizedneumannproblem_main.cc
  ${DIR}/regularizedneumannproblem.cc
  ${DIR}/getgalerkinlse.h
  ${DIR}/regularizedneumannproblem.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
