#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricfiniteelements.h
  ${DIR}/parametricfiniteelements.cc
  ${DIR}/parametricfiniteelements_main.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
