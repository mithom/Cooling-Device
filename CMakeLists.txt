cmake_minimum_required(VERSION 3.6)
project(Cooling_Device)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        heat_equation_solver/HeatSolver.cpp
        heat_equation_solver/HeatSolver.hpp
        pyplotter/pyplotter.py
        build-home.sh
        build.sh
        README.md)

add_executable(Cooling_Device ${SOURCE_FILES})