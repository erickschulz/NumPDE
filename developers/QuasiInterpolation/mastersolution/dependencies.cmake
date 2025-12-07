#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/quasiinterpolation_main.cc
  ${DIR}/quasiinterpolation.h
  ${DIR}/quasiinterpolation.cc
  ${DIR}/iohelper.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
