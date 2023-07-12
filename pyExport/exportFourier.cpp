#include "/usr/include/python3.11/Python.h" //Python.h
#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include "../src/fourier.hpp"

typedef struct
{
    PyObject_HEAD
        FourierTransform *transform;
} FourierObject;

static void Fourier_dealloc(FourierObject *self)
{
    delete self->transform;
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Fourier_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    FourierObject *self;
    self = (FourierObject *)type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->transform = NULL;
    }
    return (PyObject *)self;
}

static int Fourier_init(FourierObject *self, PyObject *args, PyObject *kwds)
{
    int nSamples;
    double deltaT;
    double t0;

    if (!PyArg_ParseTuple(args, "idd", &nSamples, &deltaT, &t0))
        return -1;

    self->transform = new FourierTransform(nSamples, deltaT, t0);

    return 0;
}

static PyObject *Fourier_forwardTransform(FourierObject *self, PyObject *args)
{
    PyObject *pyInput;

    if (!PyArg_ParseTuple(args, "O", &pyInput))
        return NULL;

    // Convert Python input to std::vector<std::complex<double>>
    std::vector<std::complex<double>> input;
    Py_ssize_t inputSize = PyList_Size(pyInput);
    for (Py_ssize_t i = 0; i < inputSize; ++i)
    {
        PyObject *pyItem = PyList_GetItem(pyInput, i);
        double real = PyComplex_RealAsDouble(pyItem);
        double imag = PyComplex_ImagAsDouble(pyItem);
        input.push_back(std::complex<double>(real, imag));
    }

    // Perform forward transform
    std::vector<std::complex<double>> result = self->transform->forwardTransform(input);

    // Convert result to Python list
    PyObject *pyResult = PyList_New(result.size());
    for (size_t i = 0; i < result.size(); ++i)
    {
        PyObject *pyItem = PyComplex_FromDoubles(std::real(result[i]), std::imag(result[i]));
        PyList_SetItem(pyResult, i, pyItem);
    }

    return pyResult;
}

static PyObject *Fourier_backwardTransform(FourierObject *self, PyObject *args)
{
    PyObject *pyInput;

    if (!PyArg_ParseTuple(args, "O", &pyInput))
        return NULL;

    // Convert Python input to std::vector<std::complex<double>>
    std::vector<std::complex<double>> input;
    Py_ssize_t inputSize = PyList_Size(pyInput);
    for (Py_ssize_t i = 0; i < inputSize; ++i)
    {
        PyObject *pyItem = PyList_GetItem(pyInput, i);
        double real = PyComplex_RealAsDouble(pyItem);
        double imag = PyComplex_ImagAsDouble(pyItem);
        input.push_back(std::complex<double>(real, imag));
    }

    // Perform backward transform
    std::vector<std::complex<double>> result = self->transform->backwardTransform(input);

    // Convert result to Python list
    PyObject *pyResult = PyList_New(result.size());
    for (size_t i = 0; i < result.size(); ++i)
    {
        PyObject *pyItem = PyComplex_FromDoubles(std::real(result[i]), std::imag(result[i]));
        PyList_SetItem(pyResult, i, pyItem);
    }

    return pyResult;
}

static PyMethodDef Fourier_methods[] = {
    {"forwardTransform", (PyCFunction)Fourier_forwardTransform, METH_VARARGS, "Perform forward Fourier transform."},
    {"backwardTransform", (PyCFunction)Fourier_backwardTransform, METH_VARARGS, "Perform backward Fourier transform."},
    {NULL, NULL, 0, NULL}};

static PyTypeObject FourierType = {
    PyVarObject_HEAD_INIT(NULL, 0) "fourier.FourierTransform", // tp_name
    sizeof(FourierObject),                                     // tp_basicsize
    0,                                                         // tp_itemsize
    (destructor)Fourier_dealloc,                               // tp_dealloc
    0,                                                         // tp_print
    0,                                                         // tp_getattr
    0,                                                         // tp_setattr
    0,                                                         // tp_reserved
    0,                                                         // tp_repr
    0,                                                         // tp_as_number
    0,                                                         // tp_as_sequence
    0,                                                         // tp_as_mapping
    0,                                                         // tp_hash
    0,                                                         // tp_call
    0,                                                         // tp_str
    0,                                                         // tp_getattro
    0,                                                         // tp_setattro
    0,                                                         // tp_as_buffer
    Py_TPFLAGS_DEFAULT,                                        // tp_flags
    "Fourier object",                                          // tp_doc
    0,                                                         // tp_traverse
    0,                                                         // tp_clear
    0,                                                         // tp_richcompare
    0,                                                         // tp_weaklistoffset
    0,                                                         // tp_iter
    0,                                                         // tp_iternext
    Fourier_methods,                                           // tp_methods
    0,                                                         // tp_members
    0,                                                         // tp_getset
    0,                                                         // tp_base
    0,                                                         // tp_dict
    0,                                                         // tp_descr_get
    0,                                                         // tp_descr_set
    0,                                                         // tp_dictoffset
    (initproc)Fourier_init,                                    // tp_init
    0,                                                         // tp_alloc
    Fourier_new                                                // tp_new
};

static PyModuleDef fouriermodule = {
    PyModuleDef_HEAD_INIT,
    "fourier",
    "Module for Fourier Transform",
    -1,
    NULL, NULL, NULL, NULL, NULL};

PyMODINIT_FUNC PyInit_fourier(void)
{
    PyObject *m;

    if (PyType_Ready(&FourierType) < 0)
        return NULL;

    m = PyModule_Create(&fouriermodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&FourierType);
    PyModule_AddObject(m, "FourierTransform", (PyObject *)&FourierType);

    return m;
}
