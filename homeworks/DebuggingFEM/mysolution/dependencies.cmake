# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/debuggingfem_main.cc
  ${DIR}/locallaplaceqfe.h
  ${DIR}/locallaplaceqfe.cc
  ${DIR}/qfeinterpolator.h
  ${DIR}/qfeinterpolator.cc
  ${DIR}/qfeprovidertester.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
