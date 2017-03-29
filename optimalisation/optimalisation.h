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
    int q = 4;
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
        cout <<"in get_info"<<endl;
        // The problem N^2 variables
        n = N*N;

        // one inequality constraint
        m = 1;

        // in this example the jacobian is dense and contains N^2 nonzeros
        nnz_jac_g = N*N;

        // the hessian of the constraints contains 0 nonzero elements TODO: dit is fout!!!!!!!
        //nnz_h_lag = N*N*N*N;

        // use the C style indexing (0-based)
        index_style = TNLP::C_STYLE;
        cout <<"done getting info"<<endl;
        return true;
    }

    // returns the variable bounds
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                         Index m, Number* g_l, Number* g_u)
    {
        cout <<"getting bounds info"<<endl;
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
        g_u[0] = 0.38;

        cout <<"done getting bounds info"<<endl;
        return true;
    }

    // returns the initial point for the problem
    bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda)
    {
        cout <<"getting starting points"<<endl;
        // Here, we assume we only have starting values for x, if you code
        // your own NLP, you can provide starting values for the dual variables
        // if you wish
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        // initialize to the given starting point

        for (int i=0;i<n;i++){
            x[i] = 0.1;
        }
        cout <<"done getting starting points"<<endl;
        return true;
    }

    // returns the value of the objective function
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        assert(n == N*N);
        if(! new_x) cout << "X did not change"<<endl;
        Eigen::SparseMatrix<double>A(n,n);

        double k[N][N];
        for(int i =0;i<n;i++) k[(int)floor(i/N)][(int)(i%N)] = 0.1 + 80*pow(x[i],3.0); //0.1 + x[i]*79.9/(1+q*(1-x[i]));
        solve_heath(k,solution,A);
        //plot(solution, n);
        double obj_value_double = (Eigen::RowVectorXd::Ones(n)+Eigen::RowVectorXd::Map(solution,n)-((1/n)*Eigen::RowVectorXd::Map(solution,n)*Eigen::VectorXd::Ones(n)*Eigen::RowVectorXd::Ones(n)))*Eigen::VectorXd::Map(solution,n);
        obj_value = (Number) obj_value_double/(2*n);
        return true;
    }

    // return the gradient of the objective function grad_{x} f(x)
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        assert(n == N*N);
        if(A.nonZeros() == 0 ){
            double k[N][N];
            for(int i =0;i<n;i++) k[(int)floor(i/N)][(int)(i%N)] = 0.1 + 80*pow((x[i]),3.0);//0.1 + x[i]*79.9/(1+q*(1-x[i]));
            solve_heath(k,solution,A);
        }
        Adjoint<N> adj(A,solution);
        adj.get_jacobi_x((double*)(grad_f),q,(double*) x);
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
        cout <<"total constraint average"<<g[0]<<endl;
        return true;
    }

    // return the structure or values of the jacobian of g
    bool eval_jac_g(Index n, const Number* x, bool new_x,
                    Index m, Index nele_jac, Index* iRow, Index *jCol,
                    Number* values) {
        if(values == NULL){
            for(int i = 0; i<n;i++){
                iRow[i] = 0;
                jCol[i]= i;
            }
        }else{
            for (int i = 0; i < n; i++) {
                values[i] = 1 / (n);
            }
        }
        assert(n==nele_jac);
        return true;
    }

    //return the structure or values of the hessian of the lagrangian:
    // sigma_f * hes(f(x)) + sum_i(lambda_i*g_i(x))
    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values)
    {
        /*if(values == NULL){
            for(int i = 0; i<n;i++){
                iRow[i] = (int)floor(i/N);
                jCol[i]= (int)(i%N);
            }
        }else{
            for(int i=0;i<n*n;i++) values[i] = 0; //TODO: benaderen!
        }*/
        return false;
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
        }
        double solution[n];
        double k[N][N];
        for(int i =0;i<n;i++) k[(int)floor(i/N)][(int)(i%N)] = 0.1 + 80*pow((x[i]),3.0);//0.1 + x[i]*79.9/(1+q*(1-x[i]));
        solve_heath(k,solution,A);
        Adjoint<N> adj(A,solution);
        double grad_f[N*N];
        adj.get_jacobi_x((double*)(grad_f),q,(double*) x);
        //cout<<Eigen::VectorXd::Map(grad_f,N*N);
        //a little python experiment
        plot(solution, n, (double*)grad_f, (double*) x);
    }
};


#endif //COOLING_DEVICE_OPTIMALISATION_H
