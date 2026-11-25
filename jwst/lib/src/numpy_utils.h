#ifndef NUMPY_UTILS_H
#define NUMPY_UTILS_H

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_21_API_VERSION
#endif

#include <numpy/arrayobject.h>

static inline PyArrayObject *
create_empty_array(int nd, npy_intp *dims, int type_num)
{
    return (PyArrayObject *) PyArray_Empty(nd, dims, PyArray_DescrFromType(type_num), 0);
}

static inline PyArrayObject *
ensure_array(PyObject *obj, int *is_copy)
{
    if (PyArray_CheckExact(obj) && PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
        PyArray_TYPE((PyArrayObject *) obj) == NPY_DOUBLE) {
        *is_copy = 0;
        return (PyArrayObject *) obj;
    } else {
        *is_copy = 1;
        return (PyArrayObject *) PyArray_FromAny(
            obj, PyArray_DescrFromType(NPY_DOUBLE), 0, 0, NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST,
            NULL);
    }
}

static inline PyArrayObject *
ensure_array_int(PyObject *obj, int *is_copy)
{
    if (PyArray_CheckExact(obj) && PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
        PyArray_TYPE((PyArrayObject *) obj) == NPY_INT) {
        *is_copy = 0;
        return (PyArrayObject *) obj;
    } else {
        *is_copy = 1;
        return (PyArrayObject *) PyArray_FromAny(
            obj, PyArray_DescrFromType(NPY_INT), 0, 0, NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST,
            NULL);
    }
}

static inline int
get_array_size(PyArrayObject *arr)
{
    return (int) PyArray_Size((PyObject *) arr);
}

static inline PyArrayObject *
array_from_vector1d(npy_intp length, int typenum, void *data)
{
    PyArrayObject *arr;
    arr = (PyArrayObject *) PyArray_SimpleNewFromData(1, &length, typenum, data);
    if (arr) {
        PyArray_ENABLEFLAGS(arr, NPY_ARRAY_OWNDATA);
    }
    return arr;
}

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif // NUMPY_UTILS_H
