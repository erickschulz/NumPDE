# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/convblfmatrixprovider.h
  ${DIR}/convblfmatrixprovider.cc
  ${DIR}/convblfmatrixprovider_main.cc
  )

set(LIBRARIES Eigen3::Eigen LF_ALL)
