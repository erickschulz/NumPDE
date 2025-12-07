#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/electrostaticforce_main.cc
  ${DIR}/electrostaticforce.h
  ${DIR}/electrostaticforce.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
