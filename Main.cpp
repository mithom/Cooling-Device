#include <iostream>
#include <IpIpoptApplication.hpp>
#include "optimalisation/optimalisation.h"
#include "pyplotter/pyplotter.hpp"
using namespace std;
int n_var;

int main(int argc, char** argv){
    const int n = 50; // moet gekend zijn om snel te zijn, anders met pointers, maar trager //http://stackoverflow.com/questions/8767166/passing-a-2d-array-to-a-c-function
    /*double k[n][n];
    std::fill_n(k[0],n,0.1);
    for(int i=1;i<n;i++){
        std::fill_n(k[i],n,0.4);
    }
    double solution[n*n];
    int m = n*n;
    Eigen::SparseMatrix<double> A(m,m);
    solve_heath(k,solution, A);
    plot(solution,n*n, *k);
    cout<<"heath solved"<<endl;
    //Adjoint<n> adj = Adjoint<n>(&A,solution);
    //double dydp[n*n];
    //cout <<"adjoint initialized"<<endl;
    //adj.get_jacobi_x(dydp);
    //cout << "jacobi calculated"<<endl;*/

    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new HS071_NLP<n>();
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetNumericValue("tol", 1e-9);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return (int) status;
    }
cout <<"starting to optimise"<<endl;
    status = app->OptimizeTNLP(mynlp);

    if (status == Ipopt::Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
    }
    else {
        printf("\n\n*** The problem FAILED!\n");
    }

    return 0;
}
