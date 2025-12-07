# Build rule for problems
function(build_problem TARGET DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

  # Create OBJECT library - compiles sources once, used by both executable and tests
  add_library(${TARGET}.obj OBJECT ${SOURCES})
  target_compile_definitions(${TARGET}.obj PRIVATE CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/${DIR}")
  target_compile_definitions(${TARGET}.obj PRIVATE CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  # Object library needs library dependencies for compilation (headers, compile flags)
  target_link_libraries(${TARGET}.obj PRIVATE ${LIBRARIES})

  # Create executable using the object library's compiled files
  add_executable(${TARGET} $<TARGET_OBJECTS:${TARGET}.obj>)
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  # Link libraries to executable (not duplicated since object library uses PRIVATE)
  target_link_libraries(${TARGET} PUBLIC ${LIBRARIES})

  # Keep .static as alias to .obj for backwards compatibility
  add_library(${TARGET}.static ALIAS ${TARGET}.obj)
endfunction(build_problem)

# Build rule for tests
function(build_test TARGET TARGET_TO_TEST DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/test/dependencies.cmake)
  include(GoogleTest)

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_compile_definitions(${TARGET} PRIVATE CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/test")
  target_compile_definitions(${TARGET} PRIVATE CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  # Add object files from main problem (compiled once, reused here)
  target_sources(${TARGET} PRIVATE $<TARGET_OBJECTS:${TARGET_TO_TEST}.obj>)
  # Link test libraries and main problem libraries (from parent dependencies.cmake)
  get_target_property(MAIN_LIBS ${TARGET_TO_TEST}.obj LINK_LIBRARIES)
  # Only link test-specific libraries, main libs already come from object files
  # Filter out libraries that are already in MAIN_LIBS to avoid duplication
  set(TEST_ONLY_LIBS ${LIBRARIES})
  foreach(lib ${MAIN_LIBS})
    list(REMOVE_ITEM TEST_ONLY_LIBS ${lib})
  endforeach()
  target_link_libraries(${TARGET} PUBLIC ${MAIN_LIBS} ${TEST_ONLY_LIBS})

  # gtest_discover_tests(${TARGET}) Not necessary given that the CI pipeline runs the tests
endfunction(build_test)

# Helper function to create relative symbolic links from the current binary directory to the source directory
function(create_relative_symlink_from_bin_dir target link_name)
  # compute relative path from current binary directory to target
  file(RELATIVE_PATH target_rel ${CMAKE_CURRENT_BINARY_DIR} ${target})
  # create symbolic links
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_rel} ${CMAKE_CURRENT_BINARY_DIR}/${link_name})
endfunction()
