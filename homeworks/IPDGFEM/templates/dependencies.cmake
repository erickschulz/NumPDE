# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/ipdgfem_main.cc
  ${DIR}/ipdgfem.h
  ${DIR}/ipdgfem.cc
)

set(LIBRARIES Eigen3::Eigen LF_ALL)
