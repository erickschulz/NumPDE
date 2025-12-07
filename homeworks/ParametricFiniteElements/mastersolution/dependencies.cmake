# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricfiniteelements.h
  ${DIR}/parametricfiniteelements.cc
  ${DIR}/parametricfiniteelements_main.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
