# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/crl2error.h
  ${DIR}/crfespace.h
  ${DIR}/nonconformingcrouzeixraviartfiniteelements.cc
  ${DIR}/nonconformingcrouzeixraviartfiniteelements.h
  ${DIR}/crl2errordirichletbvp.h
  ${DIR}/nonconformingcrouzeixraviartfiniteelements_main.cc
  ${DIR}/crdirichletbvp.h
  ${DIR}/crneumannbvp.h
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
