//
// Created by thomas on 27.03.17.
//

#ifndef COOLING_DEVICE_OPTIMALISATION_H
#define COOLING_DEVICE_OPTIMALISATION_H


// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
#include <cassert>
#include <iostream>
#include "../heath_equation_solver/HeathSolver.hpp"
#include "adjoint.h"
#include "../pyplotter/pyplotter.hpp"

using namespace std;
using namespace Ipopt;

template <int N>
class HS071_NLP : public TNLP{
    Eigen::SparseMatrix<double> A;
    double solution[N*N];
public:
    // constructor
    HS071_NLP(){
        A = Eigen::SparseMatrix<double>(N*N,N*N);
    };

    //destructor
    ~HS071_NLP(){};

    // returns the size of the problem
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                      Index& nnz_h_lag, IndexStyleEnum& index_style)
    {
        // The problem N^2 variables
        n = N*N;

        // one inequality constraint
        m = 1;

        // in this example the jacobian is dense and contains N^2 nonzeros
        nnz_jac_g = N*N;

        // the hessian of the constraints contains 0 nonzero elements TODO: dit is fout!!!!!!!
        nnz_h_lag = 0;

        // use the C style indexing (0-based)
        index_style = TNLP::C_STYLE;

        return true;
    }

    // returns the variable bounds
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                         Index m, Number* g_l, Number* g_u)
    {
        // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
        // If desired, we could assert to make sure they are what we think they are.
        assert(n == N*N);
        assert(m == 1);
        // the variables have lower bounds of 0
        for (Index i=0; i<n; i++) {
            x_l[i] = 0.0;
        }

        // the variables have upper bounds of 1
        for (Index i=0; i<n; i++) {
            x_u[i] = 1.0;
        }

        // the first constraint g1 has a lower bound of 0
        g_l[0] = 0;
        // the first constraint g1 has an upper bound of 0.4
        g_u[0] = 0.4;


        return true;
    }

    // returns the initial point for the problem
    bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda)
    {
        // Here, we assume we only have starting values for x, if you code
        // your own NLP, you can provide starting values for the dual variables
        // if you wish
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        // initialize to the given starting point

        for (int i=0;i<n;i++){
            x[i] = 0.4;
        }

        return true;
    }

    // returns the value of the objective function
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        assert(n == N*N);
        if(! new_x) cout << "X did not change"<<endl;
        Eigen::SparseMatrix<double>A(n,n);

        double k[N][N];
        for(int i =0;i<n;i++) k[(int)floor(i/N)][(int)(i%N)] = x[i]*80 + (1 - x[i])*0.01;
        solve_heath(k,solution,&A);

        for (int i=0;i<N*N;i++){
            assert(solution[i]>=300);
            obj_value += solution[i]-300;
        }

        obj_value = obj_value/(n);

        return true;
    }

    // return the gradient of the objective function grad_{x} f(x)
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        assert(n == N*N);
        Adjoint<N> adj = Adjoint<N>(&A,solution);
        double dydp[n];
        adj.get_jacobi_x((double*)(grad_f));
        return true;
    }

    // return the value of the constraints: g(x)
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
    {
        assert(n == N*N);
        assert(m == 1);

        for (int i=0;i<n;i++){
            g[0] += x[i];
        }
        g[0] = g[0]/(n);

        return true;
    }

    // return the structure or values of the jacobian of g
    bool eval_jac_g(Index n, const Number* x, bool new_x,
                    Index m, Index nele_jac, Index* iRow, Index *jCol,
                    Number* values) {
        assert(*iRow == 1);
        assert(*jCol == n);
        for (int i = 0; i < n; i++) {
            values[i] = 1 / (n);

            return true;
        }
    }

    //return the structure or values of the hessian of the lagrangian:
    // sigma_f * hes(f(x)) + sum_i(lambda_i*g_i(x))
    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values)
    {
        for(int i=0;i<n;i++) values[i] = 0; //TODO: benaderen!
        return true;
    }

    void finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq) {
        // here is where we would store the solution to variables, or write to a file, etc
        // so we could use the solution.

        // For this example, we write the solution to the console
        /*std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
        for (Index i=0; i<n; i++) {
           std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
        for (Index i=0; i<n; i++) {
          std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
        }
        for (Index i=0; i<n; i++) {
          std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Objective value" << std::endl;
        std::cout << "f(x*) = " << obj_value << std::endl;

        std::cout << std::endl << "Final value of the constraints:" << std::endl;
        for (Index i=0; i<m ;i++) {
          std::cout << "g(" << i << ") = " << g[i] << std::endl;
        }*/
        if (status != SUCCESS) {
            cout << endl << "ended with non-succes status: " << status << endl;
        } else {
            double solution[n];
            double k[N][N];
            for(int i =0;i<n;i++) k[(int)floor(i/N)][(int)(i%N)] = x[i]*80 + (1 - x[i])*0.01;
            solve_heath(k,solution,&A);
            //a little python experiment
            plot(solution, n * n);
        }
    }
};


#endif //COOLING_DEVICE_OPTIMALISATION_H