# Dependencies of mastersolution:

set(SOURCES
  ${DIR}/local_assembly.h
  ${DIR}/transpsemilagr.h
  ${DIR}/transpsemilagr.cc
  ${DIR}/transpsemilagr_main.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
