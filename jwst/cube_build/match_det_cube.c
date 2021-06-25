/*
The detector pixels are represented by a 'point could' on the sky. The IFU cube is
represented by a 3-D regular grid. This module finds the point cloud members contained
in a region centered on the center of the cube spaxel. The size of the spaxel is spatial
coordinates is cdetl1 and cdelt2, while the wavelength size is zcdelt3.
This module uses the e modified shephard weighting method to determine how to  weight each point clold member
in the spaxel.
 
Main function for Python: point_emsm

Python signature: point_emsm(nplane, cdelt1, cdelt2, zcdel3, ....)

provide more details

The output of this function is a tuple of three arrays: 
example output 

Parameters
----------
xcenters : numpy.ndarray
ycenters : numpy.ndarray
zcoord : numpy.ndarray
flux : numpy.ndarray
err : numpy.ndarray
coord1 : numpy.ndarray
coord2 : numpy.ndarray
wave : numpy.ndarray
rois_pixel : numpy.ndarray
roiw_pixel : numpy.ndarray
scalerad_pixel : numpy.ndarray

nplane : int
nwave : int
ncube : int (this is a huge number does it need to int64?)
npt : int  

cdelt1 : double
cdelt2 : double
zcdelt3: double

Returns
-------

spaxel_flux : numpy.ndarray
spaxel_weight : numpy.ndarray
spaxel_iflux : numpy.ndarray
spaxel_var : numpy.ndarray

*/

#include <stdlib.h>
#include <math.h>

#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_match_det_cube_numpy_api    //WHAT IS THIS AND WHERE IS IT USED???
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

typedef double dbl_type;
static const int npy_dbl = NPY_DOUBLE;

// if float32 accuracy is enough, use the following types instead:
/*
typedef float dbl_type;
static const int npy_dbl = NPY_FLOAT;
*/

int mem_alloc(int nelem, int **fluxv, int **weightv,  int **ifluxv) {
    int  *i;
    dbl_type *f, *w
    const char *msg = "Couldn't allocate memory for output arrays.";

    // flux:
    f = (dbl_type*) realloc(*fluxv, nelem * sizeof(dbl_type));
    if (f) {
        *fluxv = f;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // weigth:
    w = (dbl_type*) realloc(*weigthv, nelem * sizeof(dbl_type));
    if (w) {
        *weightv = w;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }
    
    // iweight
    i = (int*) realloc(*ifluxv, nelem * sizeof(int));
    if (i) {
        *ifluxv = i;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    return 0;
}

/*
Caller is responsible for de-allocating memory for the   (ASK MIHIA about this statement) 
return values: spaxel_flux, spaxel_weight, spaxel_iflux, spaxel_var
*/
int match_point_emsm(dbl_type *xc, dbl_type *yc, dbl_type *zc,
		     dbl_type *coord1, dbl_type *coord2, dbl_type *wave,
		     dbl_type *flux, dbl_type *err,
		     dbl_type *rois_pixel, dbl_type *roiw_pixel, dbl_type *scalerad_pixel,
		     int nplane, int nwave, int npt, int ncube,
		     double cdelt1, double cdelt2, double zcdelt3,
		     dbl_type **spaxel_flux, dbl_type **spaxel_weight, int **spaxel_iflux,
		     dbl_type **spaxel_var{

    int nalloc;  // number of allocated output values (not bytes)
    dbl_type *fluxv = NULL, *weightv=NULL, *varv=NULL ;  // vector for spaxel 
    int *ifluxv = NULL;  // vectors for spaxel

    // allocate memory to hold output 
    nalloc = ncube;
    if (mem_alloc(ncube, &fluxv, &weightv, &ifluxv, &varv)) return 1;

    // loop over each point cloud member and find which roi spaxels it is found
    for (k = 0; k < npt; k++) {
      
      // set up the vlaues for fluxv, weightv, ifluxv, varv
    }

    // assign output values:
    *spaxel_flux = fluxv;
    *spaxel_weight = weightv;
    *spaxel_iflux = ifluxv;
    *spaxel_var = varv;

    return 0;
}


PyArrayObject * ensure_array(PyObject *obj, int *is_copy) {
    if (PyArray_CheckExact(obj) &&
        PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
        PyArray_TYPE((PyArrayObject *) obj) == npy_dbl) {
        *is_copy = 0;
        return (PyArrayObject *) obj;
    } else {
        *is_copy = 1;
        return (PyArrayObject *) PyArray_FromAny(
            obj, PyArray_DescrFromType(npy_dbl), 0, 0,
            NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST, NULL
        );
    }
}


static PyObject * point_emsm(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xco, *yco, *zco, *fluxo, *erro, *coord1o, *coord2o, *waveo;
  PyObject *rios_pixelo, *roiw_pixelo, *scalerad_pixelo;
  
  double cdelt1, cdetl2, zcdel3;
  int nplane, nwave, npt,ncube;
  
  dbl_type *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  int *spaxel_iflux=NULL;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_rios_pixel=0, free_roiw_pixel=0, free_scalerad_pixel=0 free_flux=0, free_err=0;
  
  PyArrayObject *xc, *yc, zc, *flux, *err, *coord1, *coord2, *wave, *rios_pixel, *roiw_pixel, *scalerad_pixel;
  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_iflux_arr=NULL, *spaxel_var_arr=NULL; 
  npy_intp npy_ncube = 0;

  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOiiiiddd:point_emsm",
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro,
			&rios_pixelo, &roiw_pixelo, &scalerad_pixelo,
			&nplane, &nwave, &npt, &ncube, &cdelt1, &cdelt2, &zcdelt3)) {
    return NULL;
  }

    // check that input parameters are valid:
    if (nplane < 0) {
        PyErr_SetString(PyExc_ValueError,
            "'nplane' must be a strictly positive number.");
        return NULL;
    }

    if (nwave < 0) {
        PyErr_SetString(PyExc_ValueError,
            "'nwave' must be a strictly positive number.");
        return NULL;
    }

    if (npt < 0) {
        PyErr_SetString(PyExc_ValueError,
            "'npt' must be a strictly positive number.");
        return NULL;
    }

    if ((cdelt1 < 1) || (cdelt2 < 1)) {
        PyErr_SetString(PyExc_ValueError,
            "'cdelt1' and 'cdelt2' must be a strictly positive integer numbers.");
        return NULL;
    }
    if (zcdelt3 < 1) {
        PyErr_SetString(PyExc_ValueError,
            "'zcdelt3' must be a strictly positive integer numbers.");
        return NULL;
    }

    // ensure we are working with numpy arrays and avoid creating new ones
    // if possible:
    if ((!(xc = ensure_array(xco, &free_xc))) ||
        (!(yc = ensure_array(yco, &free_yc))) ||
	(!(zc = ensure_array(zco, &free_zc))) ||
	(!(coord1 = ensure_array(coord1o, &free_coord1))) ||
	(!(coord2 = ensure_array(coord2o, &free_coord2))) ||
	(!(wave = ensure_array(waveo, &free_wave))) ||
	(!(flux = ensure_array(fluxo, &free_flux))) ||
	(!(err = ensure_array(erroro, &free_err))) ||
	(!(rios_pixel = ensure_array(rios_pixelo, &rios_pixel))) ||
	(!(roiw_pixel = ensure_array(riow_pixelo, &roiw_pixel))) ||
	(!(scalerad_pixel = ensure_array(scalerad_pixelo, &scalerad_pixel))) ||)
      {
        goto cleanup;
    }

    n = (int) PyArray_Size((PyObject *) coord1);
    ny = (int) PyArray_Size((PyObject *) coord2);
    if (n != ny) {
        PyErr_SetString(PyExc_ValueError,
            "Input coordinate arrays of unequal size.");
        goto cleanup;
    }

    nz = (int) PyArray_Size((PyObject *) wave);
    if (!n  || !nz) {
        // 0-length input arrays. Nothing to clip. Return 0-length arrays
        spaxel_flux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, npy_dbl, 0);
        if (!spaxel_flux_arr) goto fail;

	spaxel_weight_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, npy_dbl, 0);
        if (!spaxel_weight_arr) goto fail;

	spaxel_iflux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_INT, 0);
	if (!spaxel_iflux_arr) goto fail;

	spaxel_var_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, npy_dbl, 0);
        if (!spaxel_var_arr) goto fail;

	result = Py_BuildValue("(OOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_iflux_arr, spaxel_var_arr);
	
        goto cleanup;
    }

    status = match_point_emsm(n, (dbl_type *) PyArray_DATA(xc),
			      (dbl_type *) PyArray_DATA(yc),
			      (dbl_type *) PyArray_DATA(zc),
			      (dbl_type *) PyArray_DATA(coord1),
			      (dbl_type *) PyArray_DATA(coord2),
			      (dbl_type *) PyArray_DATA(wave),			      
			      (dbl_type *) PyArray_DATA(flux),
			      (dbl_type *) PyArray_DATA(err),
			      (dbl_type *) PyArray_DATA(rois_pixel),
			      (dbl_type *) PyArray_DATA(roiw_pixel),
			      (dbl_type *) PyArray_DATA(scalerad_pixel),
			      nplane, nwave, ncube, npt, cdelt1, cdelt2, zcdelt3,
			      &spaxel_flux, &spaxel_weight, &spaxel_iflux, &spaxel_var);

    if (status) {
        goto fail;

    } else {
        // create return tuple:
        npy_ncube = (npy_intp) ncube;
	
        spaxel_flux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_ncube, npy_dbl, spaxel_flux
        );
        if (!spaxel_flux_arr) goto fail;
        spaxel_flux = NULL;

        spaxel_weight_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_ncube, npy_dbl, spaxel_weight
        );
        if (!spaxel_weight_arr) goto fail;
        spaxel_weight = NULL;

	spaxel_var_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_ncube, npy_dbl, spaxel_var
        );
        if (!spaxel_var_arr) goto fail;
        spaxel_var = NULL;

        spaxel_iflux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_cube, NPY_INT, spaxel_iflux
        );
        if (!spaxel_iflux_arr) goto fail;
        spaxel_iflux = NULL;

        PyArray_ENABLEFLAGS(spaxel_flux_arr, NPY_ARRAY_OWNDATA);
        PyArray_ENABLEFLAGS(spaxel_weight_arr, NPY_ARRAY_OWNDATA);
        PyArray_ENABLEFLAGS(spaxel_iflux_arr, NPY_ARRAY_OWNDATA);

        result = Py_BuildValue("(OOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_iflux_arr, spaxel_var_arr);
        goto cleanup;
    }

fail:
    Py_XDECREF(spaxel_flux_arr);
    Py_XDECREF(spaxel_weight_arr);
    Py_XDECREF(spaxel_iflux_arr);
    Py_XDECREF(spaxel_var_arr);
    free(spaxel_flux);
    free(spaxel_weight);
    free(spaxel_iflux);
    free(spaxel_var);

    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_MemoryError,
            "Unable to allocate memory for output arrays.");
    }

cleanup:
    if (free_xc) Py_XDECREF(xc);
    if (free_yc) Py_XDECREF(yc);
    if (free_zc) Py_XDECREF(zc);
    if (free_coord1) Py_XDECREF(coord1);
    if (free_coord2) Py_XDECREF(coord2);
    if (free_wave) Py_XDECREF(wave);
    if (free_flux) Py_XDECREF(flux);
    if (free_err) Py_XDECREF(err);
    if (free_rois_spaxel) Py_XDECREF(rois_spaxel);
    if (free_roiw_spaxel) Py_XDECREF(roiw_spaxel);
    if (free_scalerad_pixel) Py_XDECREF(scalerad_pixel);
    return result;
}


static PyMethodDef match_det_cube_methods[] =
{
    {
        "point_emsm",  point_emsm, METH_VARARGS,
        "point_emsm(put in doc string)"
    },
    {0, 0}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "match_det_cube",             /* m_name */
    "find point cloud matches for each spaxel center",  /* m_doc */
    -1,                          /* m_size */
    match_det_cube_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_match_det_cube(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
