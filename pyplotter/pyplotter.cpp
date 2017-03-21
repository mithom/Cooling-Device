#include <Python.h>
PyObject *pName, *pModule, *pDict, *pFunc;
PyObject *pArgs, *pValue;

PyObject *makelist(double array[], size_t size) {
    PyObject *l = PyList_New(size);
    for (size_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(array[i]));
    }
    return l;
}

void plot(double solution[], size_t size){
    Py_Initialize();
    pName = PyString_FromString("pyplotter.py");
    cout <<"got name"<<endl;
    if(pName != NULL) cout << "pname not null"<<endl;
    /* Error checking of pName left out */
    cout << pName<<endl;
    pModule = PyImport_Import(pName);
    cout << pModule <<endl;
    Py_DECREF(pName);
    cout <<"got file"<<endl;
    if(pModule != NULL){
        cout <<"yeay"<<endl;
    }else cout <<"oh :'("<<endl;

    pFunc = PyObject_GetAttrString(pModule, "test");
    cout <<"got func"<<endl;
    PyObject *arglist = Py_BuildValue("(S)", makelist(solution,size));
    cout << "build arglist"<<endl;
    //PyObject *result = PyObject_CallObject(pFunc, arglist);
    //cout<<"called function"<<endl;
    Py_DECREF(arglist);
    cout<<"ready to finalize"<<endl;
    Py_Finalize();
}