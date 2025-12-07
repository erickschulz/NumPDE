#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/pointevaluationrhs_main.cc
  ${DIR}/pointevaluationrhs.h
  ${DIR}/pointevaluationrhs.cc
  ${DIR}/pointevaluationrhs_norms.h
  ${DIR}/pointevaluationrhs_norms.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
