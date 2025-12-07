# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/coupledsecondorderbvp_main.cc
  ${DIR}/coupledsecondorderbvp.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
