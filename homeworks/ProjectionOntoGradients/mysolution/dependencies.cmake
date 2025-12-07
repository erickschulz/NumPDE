# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/projectionontogradients.h
  ${DIR}/projectionontogradients.cc
  ${DIR}/projectionontogradients_main.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
