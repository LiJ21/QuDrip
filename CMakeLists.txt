cmake_minimum_required(VERSION 3.10)
project(qudrip LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# output executables into a "build" folder
# use a fixed build output directory named "build" inside the project
set(BUILD_DIR ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BUILD_DIR})

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