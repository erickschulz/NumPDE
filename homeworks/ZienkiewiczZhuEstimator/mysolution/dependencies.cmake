# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/zienkiewiczzhuestimator_main.cc
  ${DIR}/zienkiewiczzhuestimator.h
  ${DIR}/zienkiewiczzhuestimator.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
