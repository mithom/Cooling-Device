#include <iostream>
#define EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_BLAS
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

const double Z = 0.001;
const int T = 300;
const double Y_r_iso = 0.003;

const int k_steel = 80;
const double k_void = 0.1;

const double Q = 2.5;

const double X = 0.01;
const double Y = X;

void test_packages(const int n){
    VectorXd vec(n);
    VectorXd v = VectorXd::Random(100);
    vec << 1,2,3,4,5,6,7,8,9,10;
    cout << "v= " << v <<endl<< "vec= " << vec << endl;
    MatrixXd m = MatrixXd::Random(100,100);
    VectorXd v1 = m.selfadjointView<Upper>().llt().solve(v);
    cout<<v1<<endl;
}

template<int n>
void solve_heath(double (&k)[n][n]){
    //test_packages(10);
    const double dx = X/n;
    const int m = n*n;
    SparseMatrix<double> A(m,m);
    A.reserve(VectorXi::Constant(m,5)); //reserve space for 5 cols/row
    VectorXd Q = VectorXd::Ones(m);
    Q=2.5/0.0001*Q;
    double mu = 1.0/(dx*dx);
    //this probably can be done with slicing in Eigen!!!! make use of this
    for(int i=1;i < n-1;i++){//0 en n-1 zijn de randgevallen
        for(int j = 1;j < n-1;j++){ //general case
            A.insert(i*n+j,(i)*n+j) = 4*k[i][j] + k[i-1][j]+ k[i+1][j]+ k[i][j-1]+ k[i][j+1];
            A.insert(i*n+j,(i-1)*n+j) = k[i][j] + k[i-1][j];
            A.insert(i*n+j,(i+1)*n+j) =k[i][j]+ k[i+1][j];
            A.insert(i*n+j,(i)*n+j+1) =k[i][j]+ k[i][j+1];
            A.insert(i*n+j,(i)*n+j-1) =k[i][j]+ k[i][j-1];
        }
    }
    for(int j = 1;j < n-1;j++){ //top edge
        A.insert((n-1)*n+j,(n-1)*n+j) = 3*k[n-1][j] + k[n-2][j]+ k[n-1][j-1]+ k[n-1][j+1];
        A.insert((n-1)*n+j,(n-2)*n+j) = k[n-1][j] + k[n-2][j];
        A.insert((n-1)*n+j,(n-1)*n+j+1) =k[n-1][j]+ k[n-1][j+1];
        A.insert((n-1)*n+j,(n-1)*n+j-1) =k[n-1][j]+ k[n-1][j-1];
    }
    for(int j = 1;j < n-1;j++){ //bottom edge
        A.insert(j,j) = 3*k[0][j] + k[1][j]+ k[0][j-1]+ k[0][j+1];
        A.insert(j,n+j) =k[0][j]+ k[1][j];
        A.insert(j,j+1) =k[0][j]+ k[0][j+1];
        A.insert(j,j-1) =k[0][j]+ k[0][j-1];
    }
    for(int i = 1;i < n-1;i++){ //right edge
        A.insert(i*n+n-1,(i)*n+n-1) = 3*k[i][n-1] + k[i-1][n-1]+ k[i+1][n-1]+ k[i][n-2];
        A.insert(i*n+n-1,(i-1)*n+n-1) = k[i][n-1] + k[i-1][n-1];
        A.insert(i*n+n-1,(i+1)*n+n-1) = k[i][n-1]+ k[i+1][n-1];
        A.insert(i*n+n-1,(i)*n+n-2) = k[i][n-1]+ k[i][n-2];
    }
    //bottom left corner
    A.insert(0,0) = 2*k[0][0] + k[1][0]+ k[0][1];
    A.insert(0,1) = k[0][0] + k[0][1];
    A.insert(0,n) =k[0][0]+ k[1][0];
    //bottom right corner
    A.insert(n-1,n-1) = 2*k[0][n-1] + k[0][n-2]+ k[1][n-1];
    A.insert(n-1,n-2) = k[0][n-1] + k[0][n-2];
    A.insert(n-1,2*n-1) =k[0][n-1]+ k[1][n-1];
    //upper left corner
    A.insert(n*(n-1),n*(n-1)) = 2*k[n-1][0] + k[n-2][0]+ k[n-1][1];
    A.insert(n*(n-1),n*(n-1)+1) = k[n-1][0] + k[n-1][1];
    A.insert(n*(n-1),n*(n-2)) =k[n-1][0]+ k[n-2][0];
    //upper right corner
    A.insert(m-1,m-1) = 2*k[n-1][n-1] + k[n-2][n-1]+ k[n-1][n-2];
    A.insert(m-1,m-2) = k[n-1][n-1] + k[n-1][n-2];
    A.insert(m-1,m-1-n) =k[n-1][n-1]+ k[n-2][n-1];
    //compensate matrix scaling
    A=-1*A*mu/2.0;

    for(int i = 1;i<n-1;i++){ //left edge
        if(i*dx <= Y_r_iso || i*dx >= X-Y_r_iso){ //upper and lower side
            A.insert(i*n,(i)*n) = (-3*k[i][0]+ k[i-1][0]+ k[i+1][0]+ k[i][1])*mu/2;
            A.insert(i*n,(i-1)*n) = (k[i][0]+ k[i-1][0])*mu/2;
            A.insert(i*n,(i+1)*n) =(k[i][0]+ k[i+1][0])*mu/2;
            A.insert(i*n,(i)*n+1) =(k[i][0]+ k[i][1])*mu/2;
        }else{ //middle
        A.insert(i*n,(i)*n) = 1;
        Q(i*n)=300;
        }

    }

    //SparseLU<SparseMatrix<double> > solver;
    /*solver.analyzePattern(A);   // for this step the numerical values of A are not used
    cout <<"analysed" <<endl;
    solver.factorize(A);
    cout << "factorised" <<endl;*/
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(A);
    cout <<"calculated" <<endl;

    VectorXd x1 = solver.solve(Q);
    cout <<"solved"<<endl;
    cout <<x1<<endl;
    cout << solver.error();
    cout<<endl<<endl<<endl<<A;

}