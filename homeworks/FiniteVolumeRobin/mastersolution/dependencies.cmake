# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/finitevolumerobin_main.cc
  ${DIR}/finitevolumerobin.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
