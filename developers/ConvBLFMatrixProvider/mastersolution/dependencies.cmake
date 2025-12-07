#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/convblfmatrixprovider.h
  ${DIR}/convblfmatrixprovider.cc
  ${DIR}/convblfmatrixprovider_main.cc
  )

set(LIBRARIES Eigen3::Eigen LF_ALL)
