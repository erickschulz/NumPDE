# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/implrk3prey_main.cc
  ${DIR}/implrk3prey.h
  ${DIR}/implrk3prey.cc
  ${DIR}/dampnewton.h
)

set(LIBRARIES
  Eigen3::Eigen
)
