#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/finitevolumerobin_main.cc
  ${DIR}/finitevolumerobin.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
