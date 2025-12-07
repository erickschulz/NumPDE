# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/advectionfv2d_main.cc
  ${DIR}/advectionfv2d.h
  ${DIR}/advectionfv2d.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
