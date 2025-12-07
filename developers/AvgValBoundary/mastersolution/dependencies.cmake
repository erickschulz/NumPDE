#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/avgvalboundary_main.cc
  ${DIR}/avgvalboundary.cc
  ${DIR}/avgvalboundary.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
