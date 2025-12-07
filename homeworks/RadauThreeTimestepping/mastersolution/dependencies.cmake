# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/radauthreetimestepping_main.cc
  ${DIR}/radauthreetimestepping.h
  ${DIR}/radauthreetimestepping.cc
  ${DIR}/radauthreetimesteppingode.h
  ${DIR}/radauthreetimesteppingode.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
