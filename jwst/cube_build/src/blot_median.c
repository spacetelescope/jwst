/*

Main function for Python: blot_wrapper

Python signature:
            result = blot_wrapper(roi_det, blot_xsize, blot_ysize, xstart, xsize2,
                                  xcenter, ycenter,
                                  x_cube, y_cube, flux_cube)

This module is used in outlier detection. Each file is mapped to the sky to create a single type
IFU cube. Using the set of single IFUs a median IFU is determined. This median image is then
mapped (blotted) back to the detector plane. This routine determines the overlap between
the median blotted back image and the detector plane.

The output of this function is a tuple of 2 arrays:(blot_flux, blot_weight)

Parameters
----------
roi : double
   Region of interest size for matching median image to detector pixels
blot_xsize : int
   X axis detector size
blot_ysize : int
   Y axis detector size
xstart : int
   Only valid for MIRI. For NIRSpec = 0
   The left most x detector pixel number the median image is being blotted to
   Blotting occurs separately for each channel. We need to know which side
   of the detector the median image is being blotted to
xsize2 : int
   If MIRI xsize2 = x size of the detector side the medina image is being blotted to
   If NIRSpec, xsize2 = xsize
xcenter : double array
   x center pixel values of the detector values. Size xsize2.
ycenter : double array
   y center pixel values of the detector values. Size ysize
x_cube : double array
   x coordinate of median cube mapped back to detector
y_cube : double array
   y coordinate of median cube mapped back to detector
flux_cube : double array
   flux of median cube

Returns
-------
blot_flux : numpy.ndarray
  IFU spaxel cflux
blot_weight : numpy.ndarray
  IFU spaxel weight
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_blot_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

//_______________________________________________________________________
// Allocate the memory to for the blotted arrays
//_______________________________________________________________________

int alloc_blot_arrays(int nelem, double **fluxv, double **weightv) {

    const char *msg = "Couldn't allocate memory for output arrays.";

    // flux:
    if (!(*fluxv  = (double*)calloc(nelem, sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, msg);
        goto failed_mem_alloc;
    }

    //weight
    if (!(*weightv  = (double*)calloc(nelem, sizeof(double)))) {
      PyErr_SetString(PyExc_MemoryError, msg);
      goto failed_mem_alloc;
    }
    return 0;

 failed_mem_alloc:
    free(*fluxv);
    free(*weightv);

}


// Match the median cube mapped to detector space with detector pixels
// Match occurs when x distance or y distance < roi size (set to 1 pixel)
// The combined flux for the blotted image is determined using a modified
// Shepard weighting methods.

// return values: blot_flux, blot_weight

int blot_overlap(double roi, int xsize_det, int ysize_det,
		 int xstart, int xsize2, int ncube,
		 double *xcenter, double *ycenter,
		 double *x_cube, double *y_cube, double *flux_cube,
		 double **blot_flux, double **blot_weight) {

  double *fluxv, *weightv;  // vector for blot values
  int k, ix, iy;
  long index2d;
  double dx, dy, dxy, weight_distance, weighted_flux;
  int npt = xsize_det * ysize_det;

  // allocate memory to hold output
  if (alloc_blot_arrays(npt, &fluxv, &weightv)) return 1;

  for (k = 0; k < ncube; k++) {
    for (ix = 0; ix< xsize2; ix ++){
      dx = fabs(x_cube[k] - xcenter[ix]);
      if( dx <= roi){
	for ( iy = 0; iy < ysize_det; iy ++){
	  dy = fabs(y_cube[k] - ycenter[iy]);
	  if (dy <= roi){
	    dxy = sqrt(dx*dx + dy*dy);
	    weight_distance = exp(-dxy);
	    weighted_flux = weight_distance * flux_cube[k];
	    index2d =  iy * xsize_det + (ix + xstart);
	    fluxv[index2d] = fluxv[index2d] + weighted_flux;
	    weightv[index2d] = weightv[index2d] + weight_distance;
	  } //
	} //end loop over iy
      } //
    } // end loop over ix
  } // end loop over ncube

  // assign output values:

  *blot_flux = fluxv;
  *blot_weight = weightv;

  return 0;
}


//  set up the C extension

PyArrayObject * ensure_array(PyObject *obj, int *is_copy) {
    if (PyArray_CheckExact(obj) &&
        PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
        PyArray_TYPE((PyArrayObject *) obj) == NPY_DOUBLE) {
        *is_copy = 0;
        return (PyArrayObject *) obj;
    } else {
        *is_copy = 1;
        return (PyArrayObject *) PyArray_FromAny(
            obj, PyArray_DescrFromType(NPY_DOUBLE), 0, 0,
            NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST, NULL
        );
    }
}


static PyObject *blot_wrapper(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xcentero, *ycentero, *x_cubeo, *y_cubeo, *flux_cubeo;

  double roi;
  int  xstart, xsize_det, ysize_det, xsize2;
  int nx, ny,nz, ncube;
  double *blot_flux=NULL, *blot_weight=NULL;
  int free_xcenter=0, free_ycenter=0, free_xcube=0, free_ycube=0, free_flux =0,status=0;

  PyArrayObject *xcenter, *ycenter, *x_cube, *y_cube, *flux_cube;
  PyArrayObject *blot_flux_arr=NULL, *blot_weight_arr=NULL;

  npy_intp npt = 0;

  if (!PyArg_ParseTuple(args, "diiiiOOOOO:blot_wrapper",
			&roi, &xsize_det,  &ysize_det, &xstart, &xsize2, &xcentero, &ycentero,
			&x_cubeo, &y_cubeo, &flux_cubeo)){
    return NULL;
  }

  // ensure we are working with numpy arrays and avoid creating new ones
  // if possible:
  if ((!(xcenter = ensure_array(xcentero, &free_xcenter))) ||
      (!(ycenter = ensure_array(ycentero, &free_ycenter))) ||
      (!(x_cube = ensure_array(x_cubeo, &free_xcube))) ||
      (!(y_cube = ensure_array(y_cubeo, &free_ycube))) ||
      (!(flux_cube = ensure_array(flux_cubeo, &free_flux)))){
    goto cleanup;
  }

  nx = (int) PyArray_Size((PyObject *) x_cube);
  ny = (int) PyArray_Size((PyObject *) y_cube);
  nz = (int) PyArray_Size((PyObject *) flux_cube);
  if (ny != nz || nz != nx ) {
    PyErr_SetString(PyExc_ValueError,
		    "Input coordinate arrays of unequal size.");
    goto cleanup;

  }

  npt = (npy_intp) (xsize_det * ysize_det);
  if (npt ==0) {
    // 0-length input arrays. Nothing to clip. Return 0-length arrays
    blot_flux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npt, NPY_DOUBLE, 0);
    if (!blot_flux_arr) goto fail;

    blot_weight_arr = (PyArrayObject*) PyArray_EMPTY(1, &npt, NPY_DOUBLE, 0);
    if (!blot_weight_arr) goto fail;

    result = Py_BuildValue("(NN)", blot_flux_arr, blot_weight_arr);

    goto cleanup;
  }

  ncube = nx;
  status = blot_overlap(roi, xsize_det, ysize_det, xstart,xsize2, ncube,
			(double *) PyArray_DATA(xcenter),
			(double *) PyArray_DATA(ycenter),
			(double *) PyArray_DATA(x_cube),
			(double *) PyArray_DATA(y_cube),
			(double *) PyArray_DATA(flux_cube),
			&blot_flux, &blot_weight );


  if (status) {
    goto fail;

  } else {
    // create return tuple:

    blot_flux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npt, NPY_DOUBLE, blot_flux);
    if (!blot_flux_arr) goto fail;
    blot_flux = NULL;

    blot_weight_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npt, NPY_DOUBLE, blot_weight);
    if (!blot_weight_arr) goto fail;
    blot_weight = NULL;

    PyArray_ENABLEFLAGS(blot_flux_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(blot_weight_arr, NPY_ARRAY_OWNDATA);

    result = Py_BuildValue("(NN)", blot_flux_arr, blot_weight_arr);
    goto cleanup;
  }

 fail:
  Py_XDECREF(blot_flux_arr);
  Py_XDECREF(blot_weight_arr);

  free(blot_flux);
  free(blot_weight);

  if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_MemoryError,
		    "Unable to allocate memory for output arrays.");
  }

 cleanup:
  if (free_xcenter)Py_XDECREF(xcenter);
  if (free_ycenter) Py_XDECREF(ycenter);
  if (free_xcube) Py_XDECREF(x_cube);
  if (free_ycube) Py_XDECREF(y_cube);
  if (free_flux) Py_XDECREF(flux_cube);

  return result;
}


static PyMethodDef blot_methods[] =
{
    {
        "blot_wrapper",
	blot_wrapper,
	METH_VARARGS,
        "blot_wrapper()"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "blot_median",             /* m_name */
    "blot median to detector space",  /* m_doc */
    -1,                          /* m_size */
    blot_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_blot_median(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
