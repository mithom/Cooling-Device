cmake_minimum_required(VERSION 3.6)
project(Cooling_Device)

set(CMAKE_CXX_STANDARD 14)
# -I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu/python2.7 -lpython2.7
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lz -Wl,-rpath,/usr/local/lib -L/usr/local/lib -lmesh_opt -I/usr/local/include -I/opt/eigen -pthread -std=gnu++11 -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization -Wno-variadic-macros -fopenmp")

set(SOURCE_FILES
        heat_equation_solver/HeatSolver.cpp
        heat_equation_solver/HeatSolver.hpp)

add_executable(Cooling_Device ${SOURCE_FILES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(Cooling_Device ${PYTHON_LIBRARIES})