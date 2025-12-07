#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/regularizedneumannproblem_main.cc
  ${DIR}/regularizedneumannproblem.cc
  ${DIR}/getgalerkinlse.h
  ${DIR}/regularizedneumannproblem.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
