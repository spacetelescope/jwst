/*
 
Main function for Python: xart_wrapper

Python signature: 
            result = xart_wrapper(imin, imax, xsize, ysize,
                 xvec, fimg, gamma, lor_amp, g_std, g_dx, g1_amp, g2_amp)

This module is used in MRS straylight subtraction.  'Straylight' is actually caused by the
detector cross-artifact arising from internal reflections within the detector substrate; in
the MIRI imager this manifests as a clear cross-shaped pattern, but in the MIRI MRS the
dispersed spectrum smears this pattern into additional traces and effectively broadens
the along-slice profile.  We fit this with a combination of a broad lorentzian profile
and four narrow gaussian profile corresponding to the cross-artifact peaks.

The parameters of these functions have been entirely determined from test data and are
encoded within a reference file; this function uses those parameters in conjunction
with the observed detector image to model the cross-artifact contribution.

The output of this function is a single array: (xart_flux) giving the 1d detector
representation of the cross-artifact model.

Parameters
----------
imin : int
   Starting column to fit (1/2 detector at a time)
imax : int
   Ending column to fit
xsize : int
   X axis detector size
ysize : int
   Y axis detector size
xvec : double array
   Vector of X pixel values across the detector.
fimg : double array
   Detector flux values mapped to a single 1d array.
gamma : double array
   Gamma parameter of the Lorenztian for each detector row.
lor_amp : double array
   Amplitude of the Lorenztian for each detector row.
g_std : double array
   Standard deviation width of the gaussians for each detector row
   1x for inner gaussian pair, 2x for outer gaussian pair.
g_dx : double array
   Linear offset of the gaussians for each detector row
   1x for inner gaussian pair, 2x for outer gaussian pair.
g1_amp : double array
    Amplitude of the inner gaussians for each detector row.
g2_amp : double array
    Amplitude of the outer gaussians for each detector row.

Returns
-------
xart_flux : numpy.ndarray
  1d cross-artifact detector model
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
      for (k = 0; k < xsize_det; k++) {
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * lor_amp[j] * gamma[j] * gamma[j])
           / (gamma[j] * gamma[j] + (xvec[k] - i) * (xvec[k] - i));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i - g_dx[j]) * (xvec[k] - i - g_dx[j])) / (2 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g1_amp[j]
           * exp(-((xvec[k] - i + g_dx[j]) * (xvec[k] - i + g_dx[j])) / (2 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g2_amp[j]
           * exp(-((xvec[k] - i - 2*g_dx[j]) * (xvec[k] - i - 2*g_dx[j])) / (8 * g_std[j] * g_std[j])));
        fluxv[xsize_det * j + k] += (fimg[j * xsize_det + i] * g2_amp[j]
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

//  Main wrapper function interface to Python code
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
