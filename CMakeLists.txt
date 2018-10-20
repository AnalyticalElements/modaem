# CMakeLists.txt - CMake configuration file for modaem
#
# CMake is a cross-platform build tool. CMake generates platform native
# makefiles that can be used with your compiler of choice. CMake uses a
# generator concept to represent different build tooling. CMake automatically
# detects the platform it is running on and generates the appropriate makefiles
# for the platform's default compiler. Different generators can also be specified.
#
# Note: CMake requires that your platform build system and compiler are
# properly installed. Build using Visual Studio requires msbuild shell.
#
# Example Usage:
#  mkdir buildproducts
#  cd buildproducts
#  cmake ..
#  make
#
# Generic Invocation:
#  cmake -E make_directory buildprod
#  cd build
#  cmake -G GENERATOR -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=1 ..
#  cmake --build . --target SOME_TARGET --config Release

cmake_minimum_required(VERSION 3.5)
project(modaem)
enable_language(Fortran)

add_subdirectory(tools)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE sources src/*.f90)

find_package(LAPACK REQUIRED) 

include_directories(${LAPACK_INCLUDE_DIR})

add_executable(modaem ${sources})

target_link_libraries(modaem ${LAPACK_LIBRARIES})
