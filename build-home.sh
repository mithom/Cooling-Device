#!/usr/bin/env bash
`libmesh-config --cxx` -I /opt/eigen/ Main.cpp `libmesh-config --cxxflags --include --ldflags --libs` -Wall -std=c++14 -llapack -lblas -llapacke `pkg-config --cflags --libs python`