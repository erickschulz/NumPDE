# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/avgvalboundary_main.cc
  ${DIR}/avgvalboundary.cc
  ${DIR}/avgvalboundary.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
