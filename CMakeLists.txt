cmake_minimum_required(VERSION 3.10)
project(qudrip LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif()

# Keep artifacts in the selected build tree. The default `-B build` layout is
# unchanged, while debug and sanitizer trees no longer overwrite it.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# common include dirs (matching your Makefile)
include_directories(${CMAKE_SOURCE_DIR}/include)
# find Eigen3 (header-only). Adjust version or PATHS if needed.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# collect all .cpp files in programs/ and add one target per file
file(GLOB PROGRAM_SOURCES "${CMAKE_SOURCE_DIR}/programs/*.cpp")

foreach(src IN LISTS PROGRAM_SOURCES)
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_options(${name} PRIVATE -O3)
  # link Eigen (adds include dirs for targets)
  target_link_libraries(${name} PRIVATE Eigen3::Eigen)
endforeach()

# regression tests (two TUs on purpose: also checks ODR/inline hygiene)
enable_testing()
add_executable(qudrip_tests
  tests/test_core.cpp
  tests/test_second_tu.cpp)
target_include_directories(qudrip_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(qudrip_tests PRIVATE Eigen3::Eigen)
add_test(NAME core COMMAND qudrip_tests)

add_executable(qudrip_cluster_tests
  tests/test_cluster.cpp
  tests/test_cluster_main.cpp)
target_include_directories(qudrip_cluster_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(qudrip_cluster_tests PRIVATE Eigen3::Eigen)
add_test(NAME cluster_constraints COMMAND qudrip_cluster_tests)

add_executable(qudrip_operator_tests
  tests/test_operator.cpp
  tests/test_operator_main.cpp)
target_include_directories(qudrip_operator_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_definitions(qudrip_operator_tests PRIVATE QUDRIP_TESTING)
target_link_libraries(qudrip_operator_tests PRIVATE Eigen3::Eigen)
add_test(NAME operator_characterization COMMAND qudrip_operator_tests)

add_executable(qudrip_index_tests tests/test_index_immutable.cpp)
target_include_directories(qudrip_index_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(qudrip_index_tests PRIVATE Eigen3::Eigen)
add_test(NAME immutable_index COMMAND qudrip_index_tests)

option(QUDRIP_BUILD_BENCHMARKS "Build opt-in QuDrip benchmarks" OFF)
if(QUDRIP_BUILD_BENCHMARKS)
  add_executable(constraint_bench benchmarks/constraint_bench.cpp)
  target_include_directories(constraint_bench PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_options(constraint_bench PRIVATE -O3)
  target_link_libraries(constraint_bench PRIVATE Eigen3::Eigen)

  add_executable(operator_bench benchmarks/operator_bench.cpp)
  target_include_directories(operator_bench PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_options(operator_bench PRIVATE -O3)
  target_link_libraries(operator_bench PRIVATE Eigen3::Eigen)
endif()
