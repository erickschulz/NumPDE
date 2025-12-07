# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/lfppdofhandling_main.cc
  ${DIR}/lfppdofhandling.h
  ${DIR}/lfppdofhandling.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
