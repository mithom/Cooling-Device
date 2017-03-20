#include <Python.h>
#include <iostream>
#include "heath_equation_solver/HeathSolver.cpp"
#include <algorithm>

using namespace std;

PyObject *pValue;
int n_var;

int main(int argc, char** argv){
    /*if(argc < 2){
        cout << "n is defaulted to 10";
        n_var=10;
    }else{
        n_var = atoi(argv[1]);
    }*/
    const int n = 10; // moet gekend zijn om snel te zijn, anders met pointers, maar trager //http://stackoverflow.com/questions/8767166/passing-a-2d-array-to-a-c-function
    double k[n][n];
    for(int i=0;i<n;i++){
        std::fill_n(k[i],n,80.0);
    }
    solve_heath(k);

    //a little python experiment
    Py_Initialize();
    int a = 5;
    pValue = PyInt_FromLong(a);
    PyRun_SimpleString("print('hello')");
    Py_Finalize();
    return 0;
}
