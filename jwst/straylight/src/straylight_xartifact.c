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



int compute_xart(float *fvec, float *xvec, float *spaxel_flux, float *lorgamma,
	       float *loramp,
	       float *gstd, float *gdx,
	       float *g1amp, float *g2amp,
	       int xsize, int ysize, int imin, int imax, long fsize) {

 long i,j,k;

  // allocate memory to hold output
 // if (alloc_flux_arrays(fsize, &fluxv)) return 1;
  printf("Now in C code");
  printf('xvec = %f, loramp = %f',xvec[10],loramp[500]);

  for (j = 500; j < 502; j++) {
   printf("j loop value %ld\n", j);
   for (i = imin; i < imax-1; i++) {
      //printf("i loop value %ld\n", i);
      for (k = imin; k < imax-1; k++) {
        spaxel_flux[xsize * j + k] += (fvec[j * xsize + i] * loramp[j] * lorgamma[j] * lorgamma[j])
           / (lorgamma[j] * lorgamma[j] + (xvec[k] - i) * (xvec[k] - i));
        //printf("k loop value %ld, xvec %f\n", k, xvec[k]);
        //spaxel_flux[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
        //   * exp(-((xvec[k] - i - gdx[j]) * (xvec[k] - i - gdx[j])) / (2 * gstd[j] * gstd[j])));
        //spaxel_flux[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
        //   * exp(-((xvec[k] - i + gdx[j]) * (xvec[k] - i + gdx[j])) / (2 * gstd[j] * gstd[j])));
       //spaxel_flux[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
        //   * exp(-((xvec[k] - i - 2*gdx[j]) * (xvec[k] - i - 2*gdx[j])) / (8 * gstd[j] * gstd[j])));
        //spaxel_flux[xsize * j + k] += (fvec[j * xsize + i] * g1amp[j]
        //   * exp(-((xvec[k] - i + 2*gdx[j]) * (xvec[k] - i + 2*gdx[j])) / (8 * gstd[j] * gstd[j])));
      }
    }
  printf("end\n");
  }

  printf("Loop done");

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

static PyObject *xartifact_wrapper(PyObject *module, PyObject *args) {


  printf("Testing 1");

  int  xsize, ysize, imin, imax, nout;
  int status=0;

  
  PyArrayObject *fvec, *xvec, *spaxel_flux, *lorgamma, *loramp, *gstd, *gdx, *g1amp, *g2amp;
  npy_intp npy_ncube = 0;

  if (!PyArg_ParseTuple(args, "OOOiiiiOOOOOO:xartifact_wrapper",
			&fvec, &xvec, &spaxel_flux, &xsize, &ysize, &imin, &imax,
			&lorgamma, &loramp, &gstd, &gdx, &g1amp, &g2amp)) {
    return NULL;
  }

  long fsize = (long) PyArray_Size((PyObject *) fvec);
  printf('fsize = %ld',fsize);


  //______________________________________________________________________
  // Do the cross-artifact calculation
  //______________________________________________________________________
  int status1 = 0;

  status = compute_xart((float *) PyArray_DATA(fvec),
			(float *) PyArray_DATA(xvec),
			(float *) PyArray_DATA(spaxel_flux),
			(float *) PyArray_DATA(lorgamma),
			(float *) PyArray_DATA(loramp),
			(float *) PyArray_DATA(gstd),
			(float *) PyArray_DATA(gdx),
			(float *) PyArray_DATA(g1amp),
			(float *) PyArray_DATA(g2amp),
			xsize, ysize, imin, imax, fsize);

   printf("Bench2");

   Py_XDECREF(fvec);
   Py_XDECREF(xvec);
   Py_XDECREF(spaxel_flux);
   Py_XDECREF(lorgamma);
   Py_XDECREF(loramp);
   Py_XDECREF(gstd);
   Py_XDECREF(gdx);
   Py_XDECREF(g1amp);
   Py_XDECREF(g2amp);

  return 0;
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
