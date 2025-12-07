# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarylength_main.cc
  ${DIR}/boundarylength.h
  ${DIR}/boundarylength.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
