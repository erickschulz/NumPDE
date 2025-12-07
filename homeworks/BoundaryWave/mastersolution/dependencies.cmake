# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarywave_main.cc
  ${DIR}/boundarywave.h
  ${DIR}/boundarywave.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
