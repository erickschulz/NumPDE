# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/linfereactdiff_main.cc
  ${DIR}/linfereactdiff.h
  ${DIR}/linfereactdiff.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
