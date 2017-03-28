//
// Created by thomas on 22.03.17.
//

#ifndef COOLING_DEVICE_ADJOINT_H
#define COOLING_DEVICE_ADJOINT_H

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Dense>
#include <math.h>
using namespace std;

template<int n>
class Adjoint {
    Eigen::RowVectorXd d;
    Eigen::MatrixXd A_inv;
    Eigen::VectorXd x;

public:
    Adjoint(Eigen::SparseMatrix<double> *A, double *x){
        using namespace Eigen;

        d = RowVectorXd::Ones(n*n);
        SparseLU<SparseMatrix<double> > solver;
        solver.analyzePattern(*A);
        solver.compute(*A);
        SparseMatrix<double> inv = SparseMatrix<double>(n*n,n*n);
        inv.reserve(VectorXd::Constant(n*n,1));
        inv.setIdentity();
        cout <<"starting to calculate"<<endl;
        A_inv = solver.solve(inv);
        cout<<"done calculating inverse"<<endl;
        (*this).x = VectorXd::Map(x,n*n);
    }

    void get_jacobi_x(double* dydp, int q){
        using namespace Eigen;

        cout << "in jacobi"<<endl;
        SparseMatrix<double> dAdpi = SparseMatrix<double>(n*n,n*n);
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

                double factor = 79.9*(1+q)/pow((q*(1-x[i*n+j])+1),2.0);
                //RowVectorXd dydpi = d * A_inv;
                //double test2 = dydpi * (dAdpi*factor) * x;
                double test2 = d*A_inv*(dAdpi*factor)*x;
                dydp[i*n+j] = test2;
            }
        }
        cout <<"general case done"<<endl;
        for(int i=1;i<n-1;i++){//left edge
            dAdpi.setZero();

            dAdpi.insert(n*i,(i)*n) = 3;
            dAdpi.insert(i*n,(i-1)*n) = -1;
            dAdpi.insert(i*n,(i+1)*n) = -1;
            dAdpi.insert(i*n,(i)*n+1) = -1;

            dAdpi.insert((i+1)*n,(i)*n) = 1;
            dAdpi.insert((i-1)*n,(i)*n) = 1;
            dAdpi.insert(i*n+1,(i)*n) = 1;

            double factor = 79.9*(1+q)/pow((q*(1-x[i*n])+1),2.0);
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * (dAdpi*factor) * x;
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
            dAdpi.insert(i*n+n-2,(i)*n+n-1) = 1;

            double factor = 79.9*(1+q)/pow((q*(1-x[i*n+n-1])+1),2.0);
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * (dAdpi*factor) * x;
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

            double factor = 79.9*(1+q)/pow((q*(1-x[i])+1),2.0);
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * (dAdpi*factor) * x;
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

            double factor = 79.9*(1+q)/pow((q*(1-x[n*(n-1)+i])+1),2.0);
            RowVectorXd dydpi = d * A_inv;
            double test2 = dydpi * (dAdpi*factor) * x;
            dydp[i*n] = test2;
        }
        //bottom left corner
        dAdpi.setZero();
        dAdpi.insert(0,0) = 2;
        dAdpi.insert(0,1) = -1;
        dAdpi.insert(0,n) = -1;

        dAdpi.insert(n,0) = 1;
        dAdpi.insert(1,0) = 1;

        double factor = 79.9*(1+q)/pow((q*(1-x[0])+1),2.0);
        RowVectorXd dydpi = d * A_inv;
        double test2 = dydpi * (dAdpi*factor) * x;
        dydp[0] = test2;

        //bottom right corner
        dAdpi.setZero();
        dAdpi.insert(n-1,n-1) = 2;
        dAdpi.insert(n-1,n-2) = -1;
        dAdpi.insert(n-1,n+n-1) = -1;

        dAdpi.insert(n+n-1,n-1) = 1;
        dAdpi.insert(n-2,n-1) = 1;

        factor = 79.9*(1+q)/pow((q*(1-x[n-1])+1),2.0);
        dydpi = d * A_inv;
        test2 = dydpi * (dAdpi*factor) * x;
        dydp[n-1] = test2;

        //top right corner
        dAdpi.setZero();
        dAdpi.insert(n*n-1,n*n-1) = 2;
        dAdpi.insert(n*n-1,n*n-2) = -1;
        dAdpi.insert(n*n-1,(n-1)*n-1) = -1;

        dAdpi.insert(n*n-2,n*n-1) = 1;
        dAdpi.insert((n-1)*n-1,n*n-1) = 1;

        factor = 79.9*(1+q)/pow((q*(1-x[n*n-1])+1),2.0);
        dydpi = d * A_inv;
        test2 = dydpi * (dAdpi*factor) * x;
        dydp[n*n-1] = test2;

        //top left corner
        dAdpi.setZero();
        dAdpi.insert((n-1)*n,(n-1)*n) = 2;
        dAdpi.insert((n-1)*n,(n-1)*n+1) = -1;
        dAdpi.insert((n-1)*n,(n-2)*n) = -1;

        dAdpi.insert((n-1)*n+1,(n-1)*n) = 1;
        dAdpi.insert((n-2)*n,(n-1)*n) = 1;

        factor = 79.9*(1+q)/pow((q*(1-x[(n-1)*n])+1),2.0);
        dydpi = d * A_inv;
        test2 = dydpi * (dAdpi*factor) * x;
        dydp[(n-1)*n] = test2;
        cout <<"adjoint calculations done"<<endl;
    }
};


#endif //COOLING_DEVICE_ADJOINT_H
