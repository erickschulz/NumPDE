# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/gausslobattoparabolic_main.cc
  ${DIR}/gausslobattoparabolic.h
  ${DIR}/gausslobattoparabolic.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
