g++ `pkg-config --cflags --libs python` -I ../packages/Eigen -I ../packages/lapack-3.7.0/LAPACKE/include HeatSolver.cpp -Wall -std=c++14 -llapack -lblas
