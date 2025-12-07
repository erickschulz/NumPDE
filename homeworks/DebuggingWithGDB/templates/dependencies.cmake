# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/debuggingwithgdb_main.cc
  ${DIR}/debuggingwithgdb.cc
  ${DIR}/debuggingwithgdb.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
