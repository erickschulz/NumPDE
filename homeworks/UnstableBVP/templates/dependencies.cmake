# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/unstablebvp_main.cc
  ${DIR}/unstablebvp.h
  ${DIR}/unstablebvp.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
