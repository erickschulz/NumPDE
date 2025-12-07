# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/nonlinschroedingerequation_main.cc
  ${DIR}/nonlinschroedingerequation.h
  ${DIR}/nonlinschroedingerequation.cc
  ${DIR}/propagator.h
  ${DIR}/propagator.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
