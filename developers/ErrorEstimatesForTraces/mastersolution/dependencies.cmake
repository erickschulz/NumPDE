#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/errorestimatesfortraces_main.cc
  ${DIR}/teelaplrobinassembly.h
  ${DIR}/teelaplrobinassembly.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
