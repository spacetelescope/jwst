
#include <Python.h>

#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static PyObject *gl_Error;

typedef npy_intp integer_t;

/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
double **ptrvector(long n)  {
    double **v;
    v=(double **)malloc((size_t) (n*sizeof(double)));
    if (!v)   {
        printf("In **ptrvector. Allocation of memory for double array failed.");
        exit(0);  }
    return v;
}


/* ==== Create Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory is allocated!                                    */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    long i,n,m;

    n=PyArray_DIMS(arrayin)[0];
    m=PyArray_DIMS(arrayin)[1];
    c=(double **)ptrvector(n);
    a=(double *) PyArray_DATA(arrayin);  /* pointer to arrayin data as double */
    for ( i=0; i<n; i++)  {
        c[i]=a+i*m;  }
    return c;
}

/* ==== Free a double *vector (vec of pointers) ========================== */
void free_Carrayptrs(double **v)  {
    free((char*) v);
}

static PyObject *
arrxyzero(PyObject *obj, PyObject *args)
{
  /* Arguments (mostly) in the order they appear */
  PyObject *oimgxy, *orefxy;
  double searchrad;

  /* Derived values */
  PyArrayObject *imgxy = NULL;
  PyArrayObject *refxy = NULL;
  PyArrayObject *ozpmat = NULL;
  double **zpmat = NULL;

  long imgnum, refnum;
  integer_t dimensions[2];
  integer_t xind, yind;
  double dx, dy;
  long j, k;

  if (!PyArg_ParseTuple(args,"OOd:arrxyzero", &oimgxy, &orefxy, &searchrad)){
    return PyErr_Format(gl_Error, "chelp.arrxyzero: Invalid Parameters.");
  }

  imgxy = (PyArrayObject *)PyArray_ContiguousFromAny(oimgxy, NPY_FLOAT32, 2, 2);
  if (!imgxy) {
    goto _exit;
  }

  refxy = (PyArrayObject *)PyArray_ContiguousFromAny(orefxy, NPY_FLOAT32, 2, 2);
  if (!refxy) {
    goto _exit;
  }

  dimensions[0] = (integer_t)(searchrad*2) + 1;
  dimensions[1] = (integer_t)(searchrad*2) + 1;
  ozpmat = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  PyArray_FILLWBYTE(ozpmat, 0);
  if (!ozpmat) {
    goto _exit;
  }
  /* Allocate memory for return matrix */
  zpmat=pymatrix_to_Carrayptrs(ozpmat);

  imgnum = PyArray_DIMS(imgxy)[0];
  refnum = PyArray_DIMS(refxy)[0];

  /* For each entry in the input image...*/
  for (j=0; j< imgnum; j++){
    /* compute the delta relative to each source in ref image */
    for (k = 0; k < refnum; k++){
        dx = *(float *)(PyArray_DATA(imgxy) + j*PyArray_STRIDES(imgxy)[0]) - *(float *)(PyArray_DATA(refxy) + k*PyArray_STRIDES(refxy)[0]);
        dy = *(float *)(PyArray_DATA(imgxy) + j*PyArray_STRIDES(imgxy)[0]+ PyArray_STRIDES(imgxy)[1]) -
             *(float *)(PyArray_DATA(refxy) + k*PyArray_STRIDES(refxy)[0]+ PyArray_STRIDES(refxy)[1]);
        if ((fabs(dx) < searchrad) && (fabs(dy) < searchrad)) {
            xind = (integer_t)(dx+searchrad);
            yind = (integer_t)(dy+searchrad);
            zpmat[yind][xind] += 1;
        }
    }
  }

 _exit:
  Py_DECREF(imgxy);
  Py_DECREF(refxy);
  free_Carrayptrs(zpmat);

  return PyArray_Return(ozpmat);
}

static PyMethodDef chelp_methods[] =
  {
    {"arrxyzero", arrxyzero, METH_VARARGS, "arrxyzero(imgxy,refxy,searchrad,zpmat)"},
    {0, 0, 0, 0}                             /* sentinel */
  };

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "chelp",
    NULL,
    -1,
    chelp_methods,
    NULL,
    NULL,
    NULL,
    NULL
};
#endif

#if PY_MAJOR_VERSION >= 3
PyObject *PyInit_chelp(void)
#else
void initchelp(void)
#endif
{
  PyObject* m;

#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  if (m == NULL) {
    return NULL;
  }

#else
  m = Py_InitModule("chelp", chelp_methods);
  if (m == NULL)
    return;
#endif

  import_array();

#if PY_MAJOR_VERSION >= 3
  return m;
#endif
}
