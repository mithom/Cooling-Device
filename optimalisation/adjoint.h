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

double harmonic_deriv(double a, double b){
    return 2*b*b/pow(a+b,2.0);
}

template<int n>
class Adjoint {
    Eigen::VectorXd d;
    //Eigen::MatrixXd A_inv;
    Eigen::VectorXd x;
    Eigen::SparseMatrix<double> A;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

public:
    Adjoint(Eigen::SparseMatrix<double> &A,double *x1){
        using namespace Eigen;

        x = VectorXd::Map(x1,n*n);
        d = (VectorXd::Ones(n*n) + (2* x) - (VectorXd::Ones(n*n)*RowVectorXd::Ones(n*n))*x*2/(n*n))/(2*n*n);

        this->A = A;
        solver.analyzePattern(this->A);
        solver.factorize(this->A);
        /*SparseMatrix<double> inv = SparseMatrix<double>(n*n,n*n);
        inv.reserve(VectorXd::Constant(n*n,1));
        inv.setIdentity();
        cout <<"starting to calculate"<<endl;
        A_inv = solver.solve(inv);
        cout<<"done calculating inverse"<<endl;*/

    }

    void get_jacobi_x(double* dydp, int q, double *k){
        using namespace Eigen;

        double mu = 1.0/((0.01/n)*(0.01/n));
        cout << "in jacobi"<<endl;
        SparseMatrix<double> dAdki = SparseMatrix<double>(n*n,n*n);
        for(int i=1;i<n-1;i++){
            for(int j=1;j<n-1;j++){//General case
                dAdki.setZero();

                dAdki.insert(i*n+j,(i)*n+j) = harmonic_deriv(k[i*n+j],k[(i-1)*n+j]) + harmonic_deriv(k[i*n+j],k[(i+1)*n+j])
                                              +harmonic_deriv(k[i*n+j],k[(i)*n+j+1]) + harmonic_deriv(k[i*n+j],k[(i)*n+j-1]);
                dAdki.insert(i*n+j,(i-1)*n+j) = -harmonic_deriv(k[i*n+j],k[(i-1)*n+j]);
                dAdki.insert(i*n+j,(i+1)*n+j) = -harmonic_deriv(k[i*n+j],k[(i+1)*n+j]);
                dAdki.insert(i*n+j,(i)*n+j+1) = -harmonic_deriv(k[i*n+j],k[(i)*n+j+1]);
                dAdki.insert(i*n+j,(i)*n+j-1) = -harmonic_deriv(k[i*n+j],k[(i)*n+j-1]);

                dAdki.insert((i+1)*n+j,(i)*n+j) = -harmonic_deriv(k[i*n+j],k[(i-1)*n+j]);
                dAdki.insert((i-1)*n+j,(i)*n+j) = -harmonic_deriv(k[i*n+j],k[(i+1)*n+j]);
                dAdki.insert(i*n+j+1,(i)*n+j) = -harmonic_deriv(k[i*n+j],k[(i)*n+j+1]);
                dAdki.insert(i*n+j-1,(i)*n+j) = -harmonic_deriv(k[i*n+j],k[(i)*n+j-1]);

                //double dkdp = 79.9*(1+q)/pow((q*(1-x[i*n+j])+1),2.0);
                double dkdp = 3*pow((x[i]),2.0);
                //RowVectorXd dydpi = d * A_inv;
                //double test2 = dydpi * (dAdki*mu*dkdp) * x;
                double test2 = d.dot(-solver.solve((dAdki*mu*dkdp)*x));
                dydp[i*n+j] = test2;
            }
        }
        cout <<"general case done"<<endl;
        for(int i=1;i<n-1;i++){//left edge
            dAdki.setZero();

            dAdki.insert(n*i,(i)*n) = harmonic_deriv(k[i*n],k[(i-1)*n]) + harmonic_deriv(k[i*n],k[(i+1)*n])
                                      +harmonic_deriv(k[i*n],k[(i)*n+1]);
            dAdki.insert(i*n,(i-1)*n) = -harmonic_deriv(k[i*n],k[(i-1)*n]);
            dAdki.insert(i*n,(i+1)*n) = -harmonic_deriv(k[i*n],k[(i+1)*n]);
            dAdki.insert(i*n,(i)*n+1) = -harmonic_deriv(k[i*n],k[(i)*n+1]);

            dAdki.insert((i+1)*n,(i)*n) = -harmonic_deriv(k[i*n],k[(i-1)*n]);
            dAdki.insert((i-1)*n,(i)*n) = -harmonic_deriv(k[i*n],k[(i+1)*n]);
            dAdki.insert(i*n+1,(i)*n) = -harmonic_deriv(k[i*n],k[(i)*n+1]);

            //double dkdp = 79.9*(1+q)/pow((q*(1-x[i*n])+1),2.0);
            double dkdp = 3*pow(x[i*n],2.0);
            //RowVectorXd dydpi = d * A_inv;
            double test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
            dydp[i*n] = test2;
        }
        for(int i=1;i<n-1;i++){//right edge
            dAdki.setZero();
            dAdki.insert(i*n+n-1,(i)*n+n-1) = harmonic_deriv(k[i*n+n-1],k[(i-1)*n+n-1]) + harmonic_deriv(k[i*n+n-1],k[(i+1)*n+n-1])
                                              + harmonic_deriv(k[i*n+n-1],k[(i)*n+n-2]);
            dAdki.insert(i*n+n-1,(i-1)*n+n-1) = -harmonic_deriv(k[i*n+n-1],k[(i-1)*n+n-1]);
            dAdki.insert(i*n+n-1,(i+1)*n+n-1) = -harmonic_deriv(k[i*n+n-1],k[(i+1)*n+n-1]);
            dAdki.insert(i*n+n-1,(i)*n+n-2) = -harmonic_deriv(k[i*n+n-1],k[(i)*n+n-2]);

            dAdki.insert((i-1)*n+n-1,(i)*n+n-1) = -harmonic_deriv(k[i*n+n-1],k[(i-1)*n+n-1]);
            dAdki.insert((i+1)*n+n-1,(i)*n+n-1) = -harmonic_deriv(k[i*n+n-1],k[(i+1)*n+n-1]);
            dAdki.insert(i*n+n-2,(i)*n+n-1) = -harmonic_deriv(k[i*n+n-1],k[(i)*n+n-2]);

            double dkdp = 3*pow((x[i*n+n-1]),2.0);
            //double dkdp = 79.9*(1+q)/pow((q*(1-x[i*n+n-1])+1),2.0);
            //RowVectorXd dydpi = d * A_inv;
            double test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
            dydp[i*n+n-1] = test2;
        }
        for(int i=1;i<n-1;i++){//bottom edge
            dAdki.setZero();

            dAdki.insert(i,(i)) = harmonic_deriv(k[i],k[i-1]) + harmonic_deriv(k[i],k[i+1]) + harmonic_deriv(k[i],k[i+n]);
            dAdki.insert(i,(i-1)) = -harmonic_deriv(k[i],k[i-1]);
            dAdki.insert(i,(i+1)) = -harmonic_deriv(k[i],k[i+1]);
            dAdki.insert(i,i+n) = -harmonic_deriv(k[i],k[i+n]);

            dAdki.insert(i-1,(i)) = -harmonic_deriv(k[i],k[i-1]);
            dAdki.insert(i+1,(i)) = -harmonic_deriv(k[i],k[i+1]);
            dAdki.insert(i+n,(i)) = -harmonic_deriv(k[i],k[i+n]);

            double dkdp = 3*pow((x[i]),2.0);
            //double dkdp = 79.9*(1+q)/pow((q*(1-x[i])+1),2.0);
            //RowVectorXd dydpi = d * A_inv;
            double test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
            dydp[i] = test2;
        }
        for(int i=1;i<n-1;i++){//top edge
            dAdki.setZero();

            dAdki.insert(n*(n-1)+i,n*(n-1)+i) = harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i]) +
                    harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i-1]) + harmonic_deriv(k[n*(n-1)+i],k[n*(n-2)+i]);
            dAdki.insert(n*(n-1)+i,n*(n-1)+i+1) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i]);
            dAdki.insert(n*(n-1)+i,n*(n-1)+i-1) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i-1]);
            dAdki.insert(n*(n-1)+i,n*(n-2)+i) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-2)+i]);

            dAdki.insert(n*(n-1)+i+1,n*(n-1)+i) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i]);
            dAdki.insert(n*(n-1)+i-1,n*(n-1)+i) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-1)+i-1]);
            dAdki.insert(n*(n-2)+i,n*(n-1)+i) = -harmonic_deriv(k[n*(n-1)+i],k[n*(n-2)+i]);

            double dkdp = 3*pow((x[n*(n-1)+i]),2.0);
            //double dkdp = 79.9*(1+q)/pow((q*(1-x[n*(n-1)+i])+1),2.0);
            //RowVectorXd dydpi = d * A_inv;
            double test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
            dydp[i*n] = test2;
        }
        //bottom left corner
        dAdki.setZero();
        dAdki.insert(0,0) = harmonic_deriv(k[0],k[1]) + harmonic_deriv(k[0],k[n]);
        dAdki.insert(0,1) = -harmonic_deriv(k[0],k[1]);
        dAdki.insert(0,n) = -harmonic_deriv(k[0],k[n]);

        dAdki.insert(n,0) = -harmonic_deriv(k[0],k[1]);
        dAdki.insert(1,0) = -harmonic_deriv(k[0],k[n]);

        double dkdp = 3*pow((x[0]),2.0);
        //double dkdp = 79.9*(1+q)/pow((q*(1-x[0])+1),2.0);
        //RowVectorXd dydpi = d * A_inv;
        double test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
        dydp[0] = test2;

        //bottom right corner
        dAdki.setZero();
        dAdki.insert(n-1,n-1) = harmonic_deriv(k[n-1],k[n-2]) + harmonic_deriv(k[n-1],k[n+n-1]);
        dAdki.insert(n-1,n-2) = -harmonic_deriv(k[n-1],k[n-2]);
        dAdki.insert(n-1,n+n-1) = -harmonic_deriv(k[n-1],k[n+n-1]);

        dAdki.insert(n+n-1,n-1) = -harmonic_deriv(k[n-1],k[n-2]);
        dAdki.insert(n-2,n-1) = -harmonic_deriv(k[n-1],k[n+n-1]);

        dkdp = 3*pow((x[n-1]),2.0);
        //dkdp = 79.9*(1+q)/pow((q*(1-x[n-1])+1),2.0);
        //dydpi = d * A_inv;
        test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
        dydp[n-1] = test2;

        //top right corner
        dAdki.setZero();
        dAdki.insert(n*n-1,n*n-1) = harmonic_deriv(k[n*n-1],k[n*n-2]) +harmonic_deriv(k[n*n-1],k[(n-1)*n-1]);
        dAdki.insert(n*n-1,n*n-2) = -harmonic_deriv(k[n*n-1],k[n*n-2]);
        dAdki.insert(n*n-1,(n-1)*n-1) = -harmonic_deriv(k[n*n-1],k[(n-1)*n-1]);

        dAdki.insert(n*n-2,n*n-1) = -harmonic_deriv(k[n*n-1],k[n*n-2]);
        dAdki.insert((n-1)*n-1,n*n-1) = -harmonic_deriv(k[n*n-1],k[(n-1)*n-1]);

        dkdp = 3*pow((x[n*n-1]),2.0);
        //dkdp = 79.9*(1+q)/pow((q*(1-x[n*n-1])+1),2.0);
        //dydpi = d * A_inv;
        test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
        dydp[n*n-1] = test2;

        //top left corner
        dAdki.setZero();
        dAdki.insert((n-1)*n,(n-1)*n) = harmonic_deriv(k[(n-1)*n],k[(n-1)*n+1]) + harmonic_deriv(k[(n-1)*n],k[(n-2)*n]);
        dAdki.insert((n-1)*n,(n-1)*n+1) = -harmonic_deriv(k[(n-1)*n],k[(n-1)*n+1]);
        dAdki.insert((n-1)*n,(n-2)*n) = -harmonic_deriv(k[(n-1)*n],k[(n-2)*n]);

        dAdki.insert((n-1)*n+1,(n-1)*n) = -harmonic_deriv(k[(n-1)*n],k[(n-1)*n+1]);
        dAdki.insert((n-2)*n,(n-1)*n) = -harmonic_deriv(k[(n-1)*n],k[(n-2)*n]);

        dkdp = 3*pow((x[(n-1)*n]),2.0);
        //dkdp = 79.9*(1+q)/pow((q*(1-x[(n-1)*n])+1),2.0);
        //dydpi = d * A_inv;
        test2 = d.dot(solver.solve((dAdki*mu*dkdp)*x));
        dydp[(n-1)*n] = test2;
        cout <<"adjoint calculations done"<<endl;
    }
};


#endif //COOLING_DEVICE_ADJOINT_H
