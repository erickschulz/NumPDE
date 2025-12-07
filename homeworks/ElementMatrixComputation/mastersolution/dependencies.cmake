# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/elementmatrixcomputation_main.cc
  ${DIR}/mylinearfeelementmatrix.h
  ${DIR}/mylinearfeelementmatrix.cc
  ${DIR}/mylinearloadvector.h
  ${DIR}/mylinearloadvector.cc
  ${DIR}/solve.h
  ${DIR}/solve.cc
  meshes/mesh.h
  meshes/mesh.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
