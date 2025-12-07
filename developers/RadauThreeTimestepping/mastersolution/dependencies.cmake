#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/radauthreetimestepping_main.cc
  ${DIR}/radauthreetimestepping.h
  ${DIR}/radauthreetimestepping.cc
  ${DIR}/radauthreetimesteppingode.h
  ${DIR}/radauthreetimesteppingode.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
