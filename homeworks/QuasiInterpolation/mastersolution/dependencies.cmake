# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/quasiinterpolation_main.cc
  ${DIR}/quasiinterpolation.h
  ${DIR}/quasiinterpolation.cc
  ${DIR}/iohelper.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
