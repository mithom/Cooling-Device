//
// Created by thomas on 27.03.17.
//

#ifndef COOLING_DEVICE_PYPLOTTER_H
#define COOLING_DEVICE_PYPLOTTER_H

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

void plot(double solution[], size_t size, double sol[], double x[]){
    Py_Initialize();
    pName = PyString_FromString("pyplotter");
    //there is no working directory since no python module, need to add this map to systpath manualy
    PyObject* sysPath = PySys_GetObject((char*)"path");
    PyObject* curDir = PyString_FromString("/home/thomas/Documents/Cooling-Device/pyplotter/");
    PyList_Append(sysPath, curDir);
    /* Error checking of pName left out */
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);
    Py_DECREF(curDir);
    if(pModule != NULL){
        cout <<"yeay"<<endl;
    }else cout <<"oh :'("<<endl;

    pFunc = PyObject_GetAttrString(pModule, "test");
    PyObject *arglist = Py_BuildValue("(S,S,S)", makelist(sol,size),makelist(x,size), makelist(sol,size));
    PyObject_CallObject(pFunc, arglist);
    cout<<"called function"<<endl;
    Py_DECREF(arglist);
    cout<<"ready to finalize"<<endl;
    Py_Finalize();
}

#endif //COOLING_DEVICE_PYPLOTTER_H
