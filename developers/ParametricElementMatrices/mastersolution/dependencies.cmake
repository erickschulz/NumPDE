#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricelementmatrices_main.cc
  ${DIR}/anisotropicdiffusionelementmatrixprovider.h
  ${DIR}/anisotropicdiffusionelementmatrixprovider.cc
  ${DIR}/fesourceelemvecprovider.h
  ${DIR}/fesourceelemvecprovider.cc
  ${DIR}/impedanceboundaryedgematrixprovider.h
  ${DIR}/impedanceboundaryedgematrixprovider.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
