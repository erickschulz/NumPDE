#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/stationarycurrents_main.cc
  ${DIR}/stationarycurrents.h
  ${DIR}/stationarycurrents.cc
  ${DIR}/stationarycurrents_supplement.h
  ${DIR}/stationarycurrents_supplement.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
