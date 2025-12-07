#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/upwindquadrature_main.cc
  ${DIR}/upwindquadrature.h
  ${DIR}/upwindquadrature.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
