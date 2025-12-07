#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/sdirkmethodoflines_main.cc
  ${DIR}/sdirkmethodoflines.h
  ${DIR}/sdirkmethodoflines.cc
  ${DIR}/sdirkmethodoflines_ode.h
  ${DIR}/sdirkmethodoflines_ode.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
