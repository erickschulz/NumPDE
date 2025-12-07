# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/testquadraturerules_main.cc
  ${DIR}/testquadraturerules.h
  ${DIR}/testquadraturerules.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
