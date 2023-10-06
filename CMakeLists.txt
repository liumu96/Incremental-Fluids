cmake_minimum_required(VERSION 3.0)

project("incremental-fluids")

set(SRC "../loadpng/lodepng.h" "../loadpng/lodepng.cpp")

add_subdirectory(1-matrixless)
add_subdirectory(2-better-advection)
add_subdirectory(3-conjugate-gradients)
add_subdirectory(4-solid-boundaries)
add_subdirectory(5-curved-boundaries)
add_subdirectory(6-heat)