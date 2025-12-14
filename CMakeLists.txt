cmake_minimum_required(VERSION 3.10)

include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/cpp-pm/hunter/archive/v0.25.3.tar.gz"
    SHA1 "0dfbc2cb5c4cf7e83533733bdfd2125ff96680cb"
    LOCAL
)

project(NPDECODES)

message("Project source dir =" "${PROJECT_SOURCE_DIR}")
message("Project binary dir =" "${PROJECT_BINARY_DIR}")

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

#setup compiler
set(CMAKE_CXX_STANDARD 20)

#enable language support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#enable warnings
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
CHECK_CXX_COMPILER_FLAG(-Wconversion COMPILER_SUPPORTS_WCONVERSION)

add_library(warning_flags INTERFACE)
if(COMPILER_SUPPORTS_WALL)
  target_compile_options(warning_flags INTERFACE -Wall)
endif()
if(COMPILER_SUPPORTS_WEXTRA)
  target_compile_options(warning_flags INTERFACE -Wextra)
endif()
if(COMPILER_SUPPORTS_WCONVERSION)
  target_compile_options(warning_flags INTERFACE -Wconversion)
endif()

# Suppress duplicate library warnings on Apple linker
if(APPLE)
  add_link_options(-Wl,-no_warn_duplicate_libraries)
endif()

# Add switches
###############################################################################
option(HOMEWORKS "usage: cmake -DHOMEWORKS=ON/OFF .." ON)
message(STATUS "HOMEWORKS ${HOMEWORKS}")

option(DEVELOPERS "usage: cmake -DDEVELOPERS=ON/OFF .." ON)
message(STATUS "DEVELOPERS ${DEVELOPERS}")

option(LECTURECODES "usage: cmake -DLECTURECODES=ON/OFF .." ON)
message(STATUS "LECTURECODES ${LECTURECODES}")

option(MASTERSOLUTION "usage: cmake -DMASTERSOLUTION=ON/OFF .." ON)
message(STATUS "MASTERSOLUTION ${MASTERSOLUTION}")

option(MYSOLUTION "usage: cmake -DMYSOLUTION=ON/OFF .." ON)
message(STATUS "MYSOLUTION ${MYSOLUTION}")

# Get Dependencies
# (don't forget to update cmake/Config.cmake.in !
###############################################################################
# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get LehrFEMpp
hunter_add_package(lehrfempp)
find_package(lehrfempp CONFIG REQUIRED)

# Create INTERFACE library that wraps all LehrFEM++ modules
# Only link top-level modules and let lower-level dependencies come in transitively
# Dependency hierarchy: uscalfe → fe → assemble → mesh → geometry → base
#
# Note: Some duplicates are unavoidable when multiple independent modules share dependencies:
# - Both lf.refinement and lf.mesh.test_utils depend on lf.mesh.utils → mesh.utils appears 2x
# - All mesh modules (hybrid2d, test_utils, refinement) depend on lf.mesh → mesh appears 2x
# The linker automatically deduplicates these - warnings are informational only
add_library(LF_ALL INTERFACE)
target_link_libraries(LF_ALL INTERFACE
  LF::lf.uscalfe
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.fe.test_utils
  LF::lf.geometry.test_utils
  LF::lf.refinement
  LF::lf.mesh.test_utils
)

# Get Google Test
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)

# Get some boost program options if LectureCodes are enabled:
if(LECTURECODES)
  hunter_add_package(Boost COMPONENTS program_options)
  find_package(Boost CONFIG REQUIRED program_options)
endif()

# include global include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add subdirectories
if(HOMEWORKS)
  add_subdirectory(homeworks)
endif()

if(DEVELOPERS)
  add_subdirectory(developers)
endif()

if(LECTURECODES)
  add_subdirectory(lecturecodes)
endif()
