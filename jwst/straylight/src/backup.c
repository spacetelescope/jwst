/*

Main function for Python: xartifact_wrapper

Python signature: result = xartifact_wrapper(fvec, xvec, xsize, ysize, imin, imax,
                                        lorgamma, loramp, gstd, gdx, g1amp, g2amp)

The output of this function is an array:(spaxel_flux)

Parameters
----------
fvec : double array
   Array of detector fluxes
xvec : double array
   x location array for indexing purposes
xsize : int
   size of detector x axis
ysize : int
   size of detector y axis
imin : int
   x axis starting element for computation
imax: int
   x axis end element for computation
lorgamma : double array
   Lorentzian function gamma (fwhm/2)
loramp : double array
   Amplitude of the Lorentzian function
gstd : double array
   1-sigma width of the gaussian function
gdx: double array
   Linear pixel offset of the gaussian function
gamp1 : double
   Amplitude of the first gaussian function
gamp2 : double
   Amplitude of the second gaussian function

Returns
-------
spaxel_flux : numpy.ndarray
  1d straylight model
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_xartifact_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

int alloc_flux_arrays(int nelem, double **fluxv) {

    const char *msg = "Couldn't allocate memory for output arrays.";

    // flux:
    if (!(*fluxv  = (double*)calloc(nelem, sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, msg);
        goto failed_mem_alloc;
    }

    return 0;

    failed_mem_alloc:
        free(*fluxv);

    return 0;
}

int compute_xart(double *fvec, double *xvec, double *lorgamma,
	       double *loramp,
	       double *gstd, double *gdx,
	       double *g1amp, double *g2amp,
	       int xsize, int ysize, int imin, int imax, long fsize,
	       double **spaxel_flux) {


  double *fluxv;  // vector for spaxel
  //fluxv[*]=0.

  long i,j,k;

  // allocate memory to hold output
  if (alloc_flux_arrays(fsize, &fluxv)) return 1;

  for (j = 0; j < ysize; j++) {
    printf("j loop value %ld\n", j);
    for (i = imin; i < imax; i++) {
      for (k = imin; k < imax; k++) {
        fluxv[xsize * j + k] += (fvec[j * xsize + i] * loramp[j] * lorgamma[j] * lorgamma[j])
           / (lorgamma[j] * lorgamma[j] + (xvec[k] - i) * (xvec[k] - i));
        fluxv[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
           * exp(-((xvec[k] - i - gdx[j]) * (xvec[k] - i - gdx[j])) / (2 * gstd[j] * gstd[j])));
        fluxv[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
           * exp(-((xvec[k] - i + gdx[j]) * (xvec[k] - i + gdx[j])) / (2 * gstd[j] * gstd[j])));
        fluxv[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
           * exp(-((xvec[k] - i - 2*gdx[j]) * (xvec[k] - i - 2*gdx[j])) / (8 * gstd[j] * gstd[j])));
        fluxv[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
           * exp(-((xvec[k] - i + 2*gdx[j]) * (xvec[k] - i + 2*gdx[j])) / (8 * gstd[j] * gstd[j])));
      }
    }
  }

  // assign output values:
  *spaxel_flux = fluxv;

  return 0;
}

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

// My inputs are fvec,xvec,xsize,ysize,imin,imax,lorgamma,loramp,gstd,gdx,g1amp,g2amp
// My output is spaxel_flux

static PyObject *xartifact_wrapper(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *fveco, *xveco, *lorgammao, *lorampo, *gstdo, *gdxo, *g1ampo, *g2ampo;

  printf("Testing 1");

  int  xsize, ysize, imin, imax, nout;
  int status=0;
  double *spaxel_flux=NULL;

  int free_fvec=0, free_xvec=0, free_lorgamma=0, free_loramp=0, free_gstd =0 , free_gdx=0;
  int free_g1amp=0, free_g2amp=0;
  
  PyArrayObject *fvec, *xvec, *lorgamma, *loramp, *gstd, *gdx, *g1amp, *g2amp;
  PyArrayObject *spaxel_flux_arr=NULL;
  npy_intp npy_ncube = 0;

  if (!PyArg_ParseTuple(args, "OOiiiiOOOOOO:xartifact_wrapper",
			&fvec, &xvec, &xsize, &ysize, &imin, &imax,
			&lorgamma, &loramp, &gstd, &gdx, &g1amp, &g2amp)) {
    return NULL;
  }


    // ensure we are working with numpy arrays and avoid creating new ones
    // if possible:
  if ((!(fvec = ensure_array(fveco, &free_fvec))) ||
      (!(xvec = ensure_array(xveco, &free_xvec))) ||
      (!(lorgamma = ensure_array(lorgammao, &free_lorgamma))) ||
      (!(loramp = ensure_array(lorampo, &free_loramp))) ||
      (!(gstd = ensure_array(gstdo, &free_gstd))) ||
      (!(gdx = ensure_array(gdxo, &free_gdx))) ||
      (!(g1amp = ensure_array(g1ampo, &free_g1amp))) ||
      (!(g2amp = ensure_array(g2ampo, &free_g2amp))) )
    {
      goto cleanup;

    }

  long fsize = (int) PyArray_Size((PyObject *) fvec);

  if ((xsize ==0)||(ysize ==0)) {
    // 0-length input arrays. Nothing to clip. Return 0-length arrays
    spaxel_flux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_DOUBLE, 0);
    if (!spaxel_flux_arr) goto fail;

    result = Py_BuildValue("(O)", spaxel_flux_arr);

    goto cleanup;

  }

  //______________________________________________________________________
  // Do the cross-artifact calculation
  //______________________________________________________________________
  int status1 = 0;
  status = compute_xart((double *) PyArray_DATA(fvec),
			(double *) PyArray_DATA(xvec),
			(double *) PyArray_DATA(lorgamma),
			(double *) PyArray_DATA(loramp),
			(double *) PyArray_DATA(gstd),
			(double *) PyArray_DATA(gdx),
			(double *) PyArray_DATA(g1amp),
			(double *) PyArray_DATA(g2amp),
			xsize, ysize, imin, imax, fsize,
			&spaxel_flux);


  if (status || status1) {
    goto fail;

  } else {
    // create return tuple:
    npy_ncube = 1;
    
    spaxel_flux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_DOUBLE, spaxel_flux);
    if (!spaxel_flux_arr) goto fail;
    spaxel_flux = NULL;

    PyArray_ENABLEFLAGS(spaxel_flux_arr, NPY_ARRAY_OWNDATA);
    result = Py_BuildValue("(O)", spaxel_flux_arr);
	
    goto cleanup;
    
  }

 fail:
  Py_XDECREF(spaxel_flux_arr);

  free(spaxel_flux);

  if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_MemoryError,
		    "Unable to allocate memory for output arrays.");
  }

 cleanup:
  if (free_fvec)Py_XDECREF(fvec);
  if (free_xvec) Py_XDECREF(xvec);
  if (free_lorgamma) Py_XDECREF(lorgamma);
  if (free_loramp) Py_XDECREF(loramp);
  if (free_gstd) Py_XDECREF(gstd);
  if (free_gdx) Py_XDECREF(gdx);
  if (free_g1amp) Py_XDECREF(g1amp);
  if (free_g2amp) Py_XDECREF(g2amp);
  return result;
}


static PyMethodDef xartifact_methods[] =
{
    {
        "xartifact_wrapper",
	xartifact_wrapper,
	METH_VARARGS,
        "xartifact_wrapper(put in doc string)"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "straylight_xartifact",             /* m_name */
    "compute cross-artifact values",  /* m_doc */
    -1,                          /* m_size */
    xartifact_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_straylight_xartifact(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
