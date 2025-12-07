# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/stableevaluationatapoint_main.cc
  ${DIR}/stableevaluationatapoint.h
  ${DIR}/stableevaluationatapoint.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
