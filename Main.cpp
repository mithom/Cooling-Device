#include <iostream>
#include "heath_equation_solver/HeathSolver.cpp"
#include <algorithm>
#include "pyplotter/pyplotter.cpp"
#include <IpIpoptNLP.hpp>
#include "optimalisation/adjoint.h"

using namespace std;

int n_var;

int main(int argc, char** argv){
    /*if(argc < 2){
        cout << "n is defaulted to 10";
        n_var=10;
    }else{
        n_var = atoi(argv[1]);
    }*/
    const int n = 100; // moet gekend zijn om snel te zijn, anders met pointers, maar trager //http://stackoverflow.com/questions/8767166/passing-a-2d-array-to-a-c-function
    double k[n][n];
    for(int i=0;i<n;i++){
        std::fill_n(k[i],n,0.1);
    }
    double solution[n*n];
    int m = n*n;
    SparseMatrix<double> A(m,m);
    solve_heath(k,solution, &A);
    cout<<"heath solved"<<endl;
    Adjoint<n> adj = Adjoint<n>(&A,solution);
    double dydp[n*n];
    cout <<"adjoint initialized"<<endl;
    adj.get_jacobi_x(dydp);
    cout << "jacobi calculated"<<endl;

    //a little python experiment
    plot(solution,n*n);
    return 0;
}
