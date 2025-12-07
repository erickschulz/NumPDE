# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/fisherkpp.cc
  ${DIR}/fisherkpp.h
  ${DIR}/fisherkpp_main.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
