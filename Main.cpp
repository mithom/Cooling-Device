#include <Python.h>
#include <iostream>
#include "heath_equation_solver/HeathSolver.cpp"

using namespace std;

PyObject *pValue;
int n;

int main(int argc, char** argv){
    if(argc < 2){
        cout << "n is defaulted to 10";
        n=10;
    }else{
        n = atoi(argv[1]);
    }
    solve_heath(n);

    //a little python experiment
    Py_Initialize();
    int a = 5;
    pValue = PyInt_FromLong(a);
    PyRun_SimpleString("print('hello')");
    Py_Finalize();
    return 0;
}
