# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/mixedfemwave_main.cc
  ${DIR}/mixedfemwave.h
  ${DIR}/mixedfemwave.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
