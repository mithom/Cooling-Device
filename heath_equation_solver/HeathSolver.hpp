#ifndef COOLING_DEVICE_HEATH_SOLVER_HPP
#define COOLING_DEVICE_HEATH_SOLVER_HPP

#include <iostream>
#define EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_BLAS
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

using namespace std;

const double Z = 0.001;
const int T = 300;
const double Y_r_iso = 0.003;

const int k_steel = 80;
const double k_void = 0.1;

const double Q = 2.5;

const double X = 0.01;
const double Y = X;

double harmonic_average(double a, double b){
    return 2/(1/a+1/b);
}

template<int n>
void solve_heath(double (*k)[n], double *result, Eigen::SparseMatrix<double> &A){
    using namespace Eigen;
    //test_packages(10);
    const double dx = X/(n-1);
    const int m = n*n;
    A.reserve(VectorXi::Constant(m,5)); //reserve space for 5 cols/row
    A.setZero();
    VectorXd Q = VectorXd::Ones(m);
    Q=2.5/(0.01*0.01)*Q;
    double mu = 1.0/(dx*dx);
    //this probably can be done with slicing in Eigen!!!! make use of this
    for(int i=1;i < n-1;i++){//0 en n-1 zijn de randgevallen
        for(int j = 1;j < n-1;j++){ //general case
            A.insert(i*n+j,(i)*n+j) = -(harmonic_average(k[i][j], k[i-1][j])+ harmonic_average(k[i][j],k[i+1][j])
                                         + harmonic_average(k[i][j],k[i][j-1])+ harmonic_average(k[i][j],k[i][j+1]));
            A.insert(i*n+j,(i-1)*n+j) = harmonic_average(k[i][j],k[i-1][j]);
            A.insert(i*n+j,(i+1)*n+j) =harmonic_average(k[i][j],k[i+1][j]);
            A.insert(i*n+j,(i)*n+j+1) =harmonic_average(k[i][j], k[i][j+1]);
            A.insert(i*n+j,(i)*n+j-1) =harmonic_average(k[i][j], k[i][j-1]);
        }
    }
    for(int j = 1;j < n-1;j++){ //top edge
        A.insert((n-1)*n+j,(n-1)*n+j) = -(harmonic_average(k[n-1][j], k[n-2][j]) +
                harmonic_average(k[n-1][j], k[n-1][j-1])+ harmonic_average(k[n-1][j],k[n-1][j+1]));
        A.insert((n-1)*n+j,(n-2)*n+j) = harmonic_average(k[n-1][j] , k[n-2][j]);
        A.insert((n-1)*n+j,(n-1)*n+j+1) =harmonic_average(k[n-1][j], k[n-1][j+1]);
        A.insert((n-1)*n+j,(n-1)*n+j-1) =harmonic_average(k[n-1][j], k[n-1][j-1]);
    }
    for(int j = 1;j < n-1;j++){ //bottom edge
        A.insert(j,j) = -(harmonic_average(k[0][j], k[1][j])
                           + harmonic_average(k[0][j],k[0][j-1])+ harmonic_average(k[0][j],k[0][j+1]));
        A.insert(j,n+j) =harmonic_average(k[0][j], k[1][j]);
        A.insert(j,j+1) =harmonic_average(k[0][j], k[0][j+1]);
        A.insert(j,j-1) =harmonic_average(k[0][j], k[0][j-1]);
    }
    for(int i = 1;i < n-1;i++){ //right edge
        A.insert(i*n+n-1,(i)*n+n-1) = -(harmonic_average(k[i][n-1],k[i-1][n-1])+
                harmonic_average(k[i][n-1],k[i+1][n-1])+ harmonic_average(k[i][n-1],k[i][n-2]));
        A.insert(i*n+n-1,(i-1)*n+n-1) = harmonic_average(k[i][n-1] , k[i-1][n-1]);
        A.insert(i*n+n-1,(i+1)*n+n-1) = harmonic_average(k[i][n-1], k[i+1][n-1]);
        A.insert(i*n+n-1,(i)*n+n-2) = harmonic_average(k[i][n-1], k[i][n-2]);
    }
    //bottom left corner
    A.insert(0,0) = -(harmonic_average(k[0][0],k[1][0])+ harmonic_average(k[0][0],k[0][1]));
    A.insert(0,1) = harmonic_average(k[0][0] , k[0][1]);
    A.insert(0,n) = harmonic_average(k[0][0], k[1][0]);
    //bottom right corner
    A.insert(n-1,n-1) = -(harmonic_average(k[0][n-1] , k[0][n-2])+ harmonic_average(k[0][n-1],k[1][n-1]));
    A.insert(n-1,n-2) = harmonic_average(k[0][n-1] , k[0][n-2]);
    A.insert(n-1,2*n-1) =harmonic_average(k[0][n-1], k[1][n-1]);
    //upper left corner
    A.insert(n*(n-1),n*(n-1)) = -(harmonic_average(k[n-1][0] , k[n-2][0])+ harmonic_average(k[n-1][0],k[n-1][1]));
    A.insert(n*(n-1),n*(n-1)+1) = harmonic_average(k[n-1][0] , k[n-1][1]);
    A.insert(n*(n-1),n*(n-2)) =harmonic_average(k[n-1][0], k[n-2][0]);
    //upper right corner
    A.insert(m-1,m-1) = -(harmonic_average(k[n-1][n-1] , k[n-2][n-1])+harmonic_average(k[n-1][n-1], k[n-1][n-2]));
    A.insert(m-1,m-2) = harmonic_average(k[n-1][n-1] , k[n-1][n-2]);
    A.insert(m-1,m-1-n) =harmonic_average(k[n-1][n-1], k[n-2][n-1]);
    //compensate matrix scaling
    (A)=-1*(A)*mu;
    for(int i = 1;i<n-1;i++){ //left edge
        if(i*dx <= Y_r_iso || i*dx >= X-Y_r_iso){ //upper and lower side
            A.insert(i*n,(i)*n) = (harmonic_average(k[i][0], k[i-1][0])+harmonic_average(k[i][0], k[i+1][0])+harmonic_average(k[i][0], k[i][1]))*mu;
            A.insert(i*n,(i-1)*n) = -(harmonic_average(k[i][0], k[i-1][0]))*mu;
            A.insert(i*n,(i+1)*n) =-(harmonic_average(k[i][0], k[i+1][0]))*mu;
            A.insert(i*n,(i)*n+1) =-(harmonic_average(k[i][0], k[i][1]))*mu;
        }else{ //middle
            A.insert(i*n,(i)*n) = 1;
            Q(i*n)=300;
        }

    }

    SparseLU<SparseMatrix<double> > solver;
    cout <<"solved"<<endl;
    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    VectorXd x1 = solver.solve(Q);
    Map<MatrixXd>( result, x1.rows(), x1.cols() ) = x1;
}

#endif //COOLING_DEVICE_HEATH_SOLVER_HPP