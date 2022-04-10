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

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_xart_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

//_______________________________________________________________________
// Allocate the memory for the output array, nelem elements long
//_______________________________________________________________________

int alloc_xart_arrays(int nelem, double **fluxv) {

    const char *msg = "Couldn't allocate memory for output arrays.";

    // flux:
    if (!(*fluxv  = (double*)calloc(nelem, sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, msg);
        goto failed_mem_alloc;
    }

    return 0;

 failed_mem_alloc:
    free(*fluxv);
   
}


// Do the computation of the cross-artifact values.

// return values: xart_flux

int xart_model(int imin, int imax, int xsize_det, int ysize_det,
		 double *xvec,
		 double *fimg, double *gamma, double *lor_amp,
		 double *g_std, double *g_dx, double *g1_amp, double *g2_amp,
		 double **xart_flux) {

  double *fluxv;  // vector for xart values
  int i, j, k;

  // This is how big the output vector needs to be
  int npt = xsize_det * ysize_det;

  // allocate memory to hold output 
  if (alloc_xart_arrays(npt, &fluxv)) return 1;

  for (j = 0; j < ysize_det; j++) {
    for (i = imin; i < imax; i++) {
      for (k = imin; k < imax; k++) {
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * lor_amp[j] * gamma[j] * gamma[j])
           / (gamma[j] * gamma[j] + (xvec[k] - i) * (xvec[k] - i));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i - g_dx[j]) * (xvec[k] - i - g_dx[j])) / (2 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i + g_dx[j]) * (xvec[k] - i + g_dx[j])) / (2 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i - 2*g_dx[j]) * (xvec[k] - i - 2*g_dx[j])) / (8 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i + 2*g_dx[j]) * (xvec[k] - i + 2*g_dx[j])) / (8 * g_std[j] * g_std[j])));
      }
    }
  }
    
  // assign output values:
  *xart_flux = fluxv;

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


static PyObject *xart_wrapper(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xveco, *fimgo, *gammao, *lor_ampo, *g_stdo, *g_dxo, *g1_ampo, *g2_ampo;

  int imin, imax;
  int xsize_det, ysize_det;
  double *xart_flux=NULL;
  int free_xvec=0, status=0;
  int free_fimg=0, free_gamma=0, free_loramp=0, free_gstd=0, free_gdx=0, free_g1amp=0, free_g2amp=0;

  PyArrayObject *xvec, *fimg, *gamma, *lor_amp, *g_std, *g_dx, *g1_amp, *g2_amp;
  PyArrayObject *xart_flux_arr=NULL;

  npy_intp npt = 0;

  if (!PyArg_ParseTuple(args, "iiiiOOOOOOOO:xart_wrapper",
			&imin, &imax, &xsize_det,  &ysize_det,
			&xveco, &fimgo, &gammao, &lor_ampo, &g_stdo, &g_dxo, &g1_ampo, &g2_ampo)){
    return NULL;
  }

  // ensure we are working with numpy arrays and avoid creating new ones
  // if possible:
  if ((!(xvec = ensure_array(xveco, &free_xvec))) ||
      (!(fimg = ensure_array(fimgo, &free_fimg))) ||
      (!(gamma = ensure_array(gammao, &free_gamma))) ||
      (!(lor_amp = ensure_array(lor_ampo, &free_loramp))) ||
      (!(g_std = ensure_array(g_stdo, &free_gstd))) ||
      (!(g_dx = ensure_array(g_dxo, &free_gdx))) ||
      (!(g1_amp = ensure_array(g1_ampo, &free_g1amp))) ||
      (!(g2_amp = ensure_array(g2_ampo, &free_g2amp)))){
    goto cleanup;
  }

  // This initializes how long the output vector should be
  npt = (npy_intp) (xsize_det * ysize_det);
  if (npt ==0) {
    // 0-length input arrays. Nothing to clip. Return 0-length arrays
    xart_flux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npt, NPY_DOUBLE, 0);
    if (!xart_flux_arr) goto fail;

    result = Py_BuildValue("(O)", xart_flux_arr);

    goto cleanup;
  }

  status = xart_model(imin, imax, xsize_det, ysize_det,
			(double *) PyArray_DATA(xvec),
			(double *) PyArray_DATA(fimg),
			(double *) PyArray_DATA(gamma),
			(double *) PyArray_DATA(lor_amp),
			(double *) PyArray_DATA(g_std),
			(double *) PyArray_DATA(g_dx),
			(double *) PyArray_DATA(g1_amp),
			(double *) PyArray_DATA(g2_amp),
			&xart_flux );


  if (status) {
    goto fail;

  } else {
    // create return tuple:
    xart_flux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npt, NPY_DOUBLE, xart_flux);
    if (!xart_flux_arr) goto fail;
    xart_flux = NULL;

    PyArray_ENABLEFLAGS(xart_flux_arr, NPY_ARRAY_OWNDATA);

    result = Py_BuildValue("(O)", xart_flux_arr);
    goto cleanup;
  }

 fail:
  Py_XDECREF(xart_flux_arr);

  free(xart_flux);

  if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_MemoryError,
		    "Unable to allocate memory for output arrays.");
  }

 cleanup:
  if (free_xvec) Py_XDECREF(xvec);
  if (free_fimg) Py_XDECREF(fimg);
  if (free_gamma) Py_XDECREF(gamma);
  if (free_loramp) Py_XDECREF(lor_amp);
  if (free_gstd) Py_XDECREF(g_std);
  if (free_gdx) Py_XDECREF(g_dx);
  if (free_g1amp) Py_XDECREF(g1_amp);
  if (free_g2amp) Py_XDECREF(g2_amp);

  return result;
}


static PyMethodDef xart_methods[] =
{
    {
        "xart_wrapper",
	xart_wrapper,
	METH_VARARGS,
        "xart_wrapper()"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "calc_xart",             /* m_name */
    "Calc x artifact",  /* m_doc */
    -1,                          /* m_size */
    xart_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_calc_xart(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
