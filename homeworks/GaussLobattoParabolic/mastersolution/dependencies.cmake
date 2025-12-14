# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/gausslobattoparabolic.h
  ${DIR}/gausslobattoparabolic.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
