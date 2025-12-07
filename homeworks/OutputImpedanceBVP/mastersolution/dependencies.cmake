# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/outputimpedancebvp_main.cc
  ${DIR}/outputimpedancebvp.h
  ${DIR}/outputimpedancebvp.cc
  ${DIR}/evalclass.h
  ${DIR}/evalclass.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
