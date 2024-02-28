cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(moving_window)
add_subdirectory(pybind11)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(moving_window
   	moving_window.cpp)
pybind11_add_module(moving_window ${moving_window})

