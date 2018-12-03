# CMakeLists.txt - CMake configuration file for modaem/tools

cmake_minimum_required(VERSION 3.5)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB_RECURSE tr0bln_sources tr0bln.f90 ../src/u_constants.f90 ../src/u_io.f90) 

file(GLOB_RECURSE contour_sources contour.f90 contouring.f90 ../src/u_constants.f90)

link_directories(/usr/local/lib)

add_executable(contour ${contour_sources})
add_executable(tr0bln ${tr0bln_sources})

target_link_libraries(contour dislin)
