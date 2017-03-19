#include <iostream>
#define EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_BLAS
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

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

void solve_heath(const int n){
    VectorXd vec(n);
    VectorXd v = VectorXd::Random(100);
    vec << 1,2,3,4,5,6,7,8,9,10;
    cout << "v= " << v <<endl<< "vec= " << vec << endl;
    MatrixXd m = MatrixXd::Random(100,100);
    VectorXd v1 = m.selfadjointView<Upper>().llt().solve(v);
    cout<<v1<<endl;
}