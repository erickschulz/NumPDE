# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/waveabc2d.h
  ${DIR}/waveabc2d.cc
  ${DIR}/waveabc2d_main.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
