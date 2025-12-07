#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarywave_main.cc
  ${DIR}/boundarywave.h
  ${DIR}/boundarywave.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
