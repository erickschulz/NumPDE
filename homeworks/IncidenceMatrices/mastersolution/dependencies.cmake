# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/incidencematrices_main.cc
  ${DIR}/incidencematrices.h
  ${DIR}/incidencematrices.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
