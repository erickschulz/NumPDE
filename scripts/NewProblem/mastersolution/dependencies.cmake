#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# C++ sources for this problem 
# DIR will be provided by the calling file.
set(SOURCES
  ${DIR}/newproblem_main.cc
  ${DIR}/newproblem.h
  ${DIR}/newproblem.cc
)

# Libraries to be used. If the code does not rely on LehrFEM++
# all the libraries LF:* can be removed 
set(LIBRARIES Eigen3::Eigen LF_ALL)
