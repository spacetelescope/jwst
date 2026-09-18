#ifndef CUBE_UTILS_H
#define CUBE_UTILS_H

#include <Python.h>

// Ensure PY_ARRAY_UNIQUE_SYMBOL is set BEFORE numpy header includes
#ifndef PY_ARRAY_UNIQUE_SYMBOL
#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_build_numpy_api
#endif

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

// Array utility declarations
PyArrayObject *
ensure_array(PyObject *obj, int *is_copy);
PyArrayObject *
ensure_array_int(PyObject *obj, int *is_copy);

// Memory allocation function declarations
int
alloc_flux_arrays(int nelem, double **fluxv, double **weightv, double **varv, double **ifluxv);

int
alloc_flux_dq_arrays(
    int nelem, double **fluxv, double **weightv, double **varv, double **ifluxv, int **dqv);

// Geometry and overlap function declarations
double
find_area_quad(double MinX, double MinY, double Xcorner[], double Ycorner[]);

double
find_area_poly(int nVertices, double xPixel[], double yPixel[]);

// Geometry / Overlap declarations
double
sh_find_overlap(
    double xcenter, double ycenter, double xlength, double ylength, double xPixelCorner[],
    double yPixelCorner[]);

#endif /* CUBE_UTILS_H */
