#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/gausslobattoparabolic_main.cc
  ${DIR}/gausslobattoparabolic.h
  ${DIR}/gausslobattoparabolic.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
