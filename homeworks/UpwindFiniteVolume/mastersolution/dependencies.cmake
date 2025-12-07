# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/upwindfinitevolume_main.cc
  ${DIR}/upwindfinitevolume.h
  ${DIR}/upwindfinitevolume.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
