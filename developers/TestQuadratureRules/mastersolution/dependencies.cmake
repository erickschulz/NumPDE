#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/testquadraturerules_main.cc
  ${DIR}/testquadraturerules.h
  ${DIR}/testquadraturerules.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
