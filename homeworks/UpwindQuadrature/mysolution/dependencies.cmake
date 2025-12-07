# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/upwindquadrature_main.cc
  ${DIR}/upwindquadrature.h
  ${DIR}/upwindquadrature.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
