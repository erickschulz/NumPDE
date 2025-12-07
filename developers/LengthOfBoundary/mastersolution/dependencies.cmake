#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarylength_main.cc
  ${DIR}/boundarylength.h
  ${DIR}/boundarylength.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
