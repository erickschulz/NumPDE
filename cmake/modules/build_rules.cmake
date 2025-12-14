# Build rule for problems
function(build_problem TARGET DIR OUTPUT_NAME)
  # Load custom dependencies if they exist
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/dependencies.cmake)
    include(${DIR}/dependencies.cmake)
  endif()

  # Auto-discover sources if not specified
  if(NOT SOURCES)
    file(GLOB MAIN_SOURCES "${DIR}/*_main.cc")
    file(GLOB OTHER_SOURCES "${DIR}/*.cc" "${DIR}/*.h")
    list(REMOVE_ITEM OTHER_SOURCES ${MAIN_SOURCES})
    set(SOURCES ${OTHER_SOURCES})
  else()
    # If custom SOURCES provided, discover main sources separately
    file(GLOB MAIN_SOURCES "${DIR}/*_main.cc")
  endif()

  # Use default libraries if not specified
  if(NOT LIBRARIES)
    set(LIBRARIES Eigen3::Eigen LF_ALL)
  endif()

  # Create object library containing core implementation
  # Compiled once and reused by both executable and test targets
  add_library(${TARGET}.obj OBJECT ${SOURCES})
  target_compile_definitions(${TARGET}.obj PRIVATE CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/${DIR}")
  target_compile_definitions(${TARGET}.obj PRIVATE CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  target_link_libraries(${TARGET}.obj PUBLIC ${LIBRARIES})

  # Create main executable by linking object library with main file
  add_executable(${TARGET} ${MAIN_SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_compile_definitions(${TARGET} PRIVATE CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/${DIR}")
  target_compile_definitions(${TARGET} PRIVATE CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  target_link_libraries(${TARGET} PRIVATE ${TARGET}.obj)
endfunction(build_problem)

# Build rule for tests
function(build_test TARGET TARGET_TO_TEST DIR OUTPUT_NAME)
  # Load custom test dependencies if they exist
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/test/dependencies.cmake)
    include(${DIR}/test/dependencies.cmake)
  endif()
  include(GoogleTest)

  # Auto-discover test sources if not specified
  if(NOT SOURCES)
    file(GLOB SOURCES "${DIR}/test/*.cc")
  endif()

  # Use default test libraries if not specified
  if(NOT LIBRARIES)
    set(LIBRARIES GTest::gtest_main)
  endif()

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_compile_definitions(${TARGET} PRIVATE CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/test")
  target_compile_definitions(${TARGET} PRIVATE CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")

  # Link against the problem's object library (inherits all dependencies automatically)
  # and add test-specific libraries
  target_link_libraries(${TARGET} PRIVATE ${TARGET_TO_TEST}.obj ${LIBRARIES})
endfunction(build_test)

# Create relative symbolic link from binary directory to source directory
# Used to make mesh files and scripts accessible from build directory
function(create_relative_symlink_from_bin_dir target link_name)
  file(RELATIVE_PATH target_rel ${CMAKE_CURRENT_BINARY_DIR} ${target})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_rel} ${CMAKE_CURRENT_BINARY_DIR}/${link_name})
endfunction()
