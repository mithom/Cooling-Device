#!/usr/bin/env bash
`libmesh-config --cxx` `pkg-config --cflags --libs python` -I /opt/eigen/ heat_equation_solver/HeatSolver.cpp `libmesh-config --cxxflags --include --ldflags --libs` -Wall -std=c++14 -llapack -lblas -llapacke -lpython2.7