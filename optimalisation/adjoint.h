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
        for(int i=1;i<n-1;i++){//left edge
            dAdpi.setZero();

            dAdpi.insert(n*i,(i)*n) = 3;
            dAdpi.insert(i*n,(i-1)*n) = -1;
            dAdpi.insert(i*n,(i+1)*n) = -1;
            dAdpi.insert(i*n,(i)*n+1) = -1;

            dAdpi.insert((i+1)*n,(i)*n) = 1;
            dAdpi.insert((i-1)*n,(i)*n) = 1;
            dAdpi.insert(i*n+1,(i)*n) = 1;
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * dAdpi * x;
            dydp[i*n] = test2;
        }
        for(int i=1;i<n-1;i++){//right edge
            dAdpi.setZero();

            dAdpi.insert(i*n+n-1,(i)*n+n-1) = 3;
            dAdpi.insert(i*n+n-1,(i-1)*n+n-1) = -1;
            dAdpi.insert(i*n+n-1,(i+1)*n+n-1) = -1;
            dAdpi.insert(i*n+n-1,(i)*n+n-2) = -1;

            dAdpi.insert((i-1)*n+n-1,(i)*n+n-1) = 1;
            dAdpi.insert((i+1)*n+n-1,(i)*n+n-1) = 1;
            dAdpi.insert(i*n+n-1,(i)*n+n-1) = 1;
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * dAdpi * x;
            dydp[i*n+n-1] = test2;
        }
        for(int i=1;i<n-1;i++){//bottom edge
            dAdpi.setZero();

            dAdpi.insert(i,(i)) = 3;
            dAdpi.insert(i,(i-1)) = -1;
            dAdpi.insert(i,(i+1)) = -1;
            dAdpi.insert(i,i+n) = -1;

            dAdpi.insert(i-1,(i)) = 1;
            dAdpi.insert(i+1,(i)) = 1;
            dAdpi.insert(i+n,(i)) = 1;
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * dAdpi * x;
            dydp[i] = test2;
        }
        for(int i=1;i<n-1;i++){//top edge
            dAdpi.setZero();

            dAdpi.insert(n*(n-1)+i,n*(n-1)+i) = 3;
            dAdpi.insert(n*(n-1)+i,n*(n-1)+i+1) = -1;
            dAdpi.insert(n*(n-1)+i,n*(n-1)+i-1) = -1;
            dAdpi.insert(n*(n-1)+i,n*(n-2)+i) = -1;

            dAdpi.insert(n*(n-1)+i+1,n*(n-1)+i) = 1;
            dAdpi.insert(n*(n-1)+i-1,n*(n-1)+i) = 1;
            dAdpi.insert(n*(n-2)+i,n*(n-1)+i) = 1;
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * dAdpi * x;
            dydp[i*n] = test2;
        }
    }
};


#endif //COOLING_DEVICE_ADJOINT_H
