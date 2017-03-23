//
// Created by thomas on 22.03.17.
//

#ifndef COOLING_DEVICE_ADJOINT_H
#define COOLING_DEVICE_ADJOINT_H

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

template<int n>
class Adjoint {
    RowVectorXd d;
    MatrixXd A_inv;
    VectorXd x;

public:
    Adjoint(SparseMatrix<double> *A, double *x){
        d = RowVectorXd::Ones(n*n);
        SparseLU<SparseMatrix<double> > solver;
        solver.analyzePattern(*A);
        solver.compute(*A);
        SparseMatrix<double> inv = SparseMatrix<double>(n*n,n*n);
        inv.reserve(VectorXd::Constant(n*n,1));
        inv.setIdentity();
        A_inv = solver.solve(inv);
        (*this).x = VectorXd::Map(x,n*n);
    }

    void get_jacobi_x(double* dydp){
        cout << "in jacobi"<<endl;
        SparseMatrix<double> dAdpi = SparseMatrix<double>(n*n,n*n);
        cout << "sparse reserved"<<endl;
        cout << "reserved"<<endl;
        for(int i=1;i<n-1;i++){
            for(int j=1;j<n-1;j++){//General case
                dAdpi.setZero();

                dAdpi.insert(i*n+j,(i)*n+j) = 4;
                dAdpi.insert(i*n+j,(i-1)*n+j) = -1;
                dAdpi.insert(i*n+j,(i+1)*n+j) = -1;
                dAdpi.insert(i*n+j,(i)*n+j+1) = -1;
                dAdpi.insert(i*n+j,(i)*n+j-1) = -1;

                dAdpi.insert((i+1)*n+j,(i)*n+j) = 1;
                dAdpi.insert((i-1)*n+j,(i)*n+j) = 1;
                dAdpi.insert(i*n+j+1,(i)*n+j) = 1;
                dAdpi.insert(i*n+j-1,(i)*n+j) = 1;
                RowVectorXd dydpi = d * A_inv;
                double test2 = dydpi * dAdpi * x;
                //double dydpi = ((*this).d)*((*this).A_inv)*(dAdpi)*((*this).x);
                dydp[i*n+j] = test2;
            }
        }
        //TODO: edge cases
    }
};


#endif //COOLING_DEVICE_ADJOINT_H
