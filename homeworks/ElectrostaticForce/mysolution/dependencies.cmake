# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/electrostaticforce_main.cc
  ${DIR}/electrostaticforce.h
  ${DIR}/electrostaticforce.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
