/*
This module forms the IFU cube in the "detector plane". This method of creating cubes is
for an engineering mode. The detector pixels are mapped to the along slice coordinate
and wavelength. The slice number represents the across slice coordinate. The pixels from
each slice are assumed to not overlap on the output plane with pixels from other slices. This
breaks the 3-D representation of the pixel's coordinates in the output plane to 2D. Standard
drizzling type routines are used to calculate the percentage of each pixel's overlap with
the output spaxels. The percentage of the pixel area  overlapped onto an output
spaxel is used to weight the pixel flux for determining the spaxel flux.


Main function for Python: cube_wrapper_internal

Python signature:
    result = cube_wrapper_internal(instrument_no, naxis1, naxis2,
                                   crval_along, cdelt_along, crval3, cdelt3,
                                   a1, a2, a3, a4, lam1, lam2, lam3, lam4,
                                   acoord, zcoord, ss,
                                   pixel_flux, pixel_err)
provide more details

The output of this function is a tuple of 5 arrays:(spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq)
example output

Parameters
----------
instrument: int
    0 = MIRI, 1 = NIRSPEC. Used for set the dq plane
naxis1 : int
   axis 1 of IFU cube
nasix2 : int
  axis 2 of IFU cube
crval_along : double
  along slice reference value in IFU cube
cdelt_along : double
  along slice sampling in IFU cube
crval3 : double
   wavelength reference value in IFU cube
cdetl3 : double
   wavelength sampling in IFU cube
a1, a2, a3, a4: double array
    array of corners of pixels holding along slice coordinates
    a1 holds corner 1
    a2 holds corner 2
    a3 holds corner 3
    a4 holds corner 4
lam1, lam2, lam3, lam4: double array
    array of corners of pixels holding wavelength coordinates
    lam1 holds corner 1
    lam2 holds corner 2
    lam3 holds corner 3
    lam4 holds corner 4
acoord : double array
   Array holds along slice coordinates in IFU cube
zcoord : double array
   Array holds the wavelength coordinates in the IFU cube
ss : double array
   Array holds the slice number
pixel_flux : double array
  Array that holds detector pixel flux
pixel_error : double array
  Array that holds detector pixel flux error

Returns
-------

spaxel_flux : array
spaxel_weight : array
spaxel_iflux : array
spaxel_var : array
spaxel_dq : array

*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_match_internal_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


// routines used from cube_utils.c
extern double sh_find_overlap(const double xcenter, const double ycenter,
		      const double xlength, const double ylength,
			      double xPixelCorner[],double yPixelCorner[]);


extern double find_area_quad(double MinX, double MinY, double Xcorner[], double Ycorner[]);


extern int alloc_flux_arrays(int nelem, double **fluxv, double **weightv, double **varv,  double **ifluxv);


//_______________________________________________________________________
//  based on a 2-D drizzling type of algorithm find the contribution of
// each detector flux to the IFU cube
//_______________________________________________________________________

int match_detector_cube(int instrument, int naxis1, int naxis2, int nz, int npt, int ncube, int na,
			double crval_along, double cdelt_along, double crval3, double cdelt3,
			double *a1, double *a2,double *a3, double*a4,
			double *lam1, double *lam2, double*lam3, double *lam4,
			double *acoord, double *zcoord, int ss,
			double *pixel_flux, double *pixel_err,
			double **spaxel_flux, double **spaxel_weight, double **spaxel_var,double **spaxel_iflux){

  double *fluxv=NULL, *weightv=NULL, *varv=NULL, *ifluxv=NULL;  // vectors for spaxel
  double along_corner[4], wave_corner[4], along_min, wave_min, along_max, wave_max, Area, MinW, MaxW, zcenter, acenter, area_overlap,
    AreaRatio, err;
  int ipixel, ia1, ia2, iz1, iz2, nplane, zz, istart, aa, j, cube_index;

  // allocate memory to hold output
  if (alloc_flux_arrays(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;


  // loop over each valid point on detector and find match to IFU cube based
  // on along slice coordinate and wavelength
  for (ipixel= 0; ipixel< npt; ipixel++){

    along_corner[0] = a1[ipixel];
    along_corner[1] = a2[ipixel];
    along_corner[2] = a3[ipixel];
    along_corner[3] = a4[ipixel];

    wave_corner[0] = lam1[ipixel];
    wave_corner[1] = lam2[ipixel];
    wave_corner[2] = lam3[ipixel];
    wave_corner[3] = lam4[ipixel];

    along_min = along_corner[0];
    wave_min = wave_corner[0];
    along_max = along_min;
    wave_max = wave_min;
    for (j = 1; j< 4; j++){
      if(along_corner[j] < along_min) { along_min = along_corner[j];}
      if(along_corner[j] > along_max) { along_max = along_corner[j];}
      if(wave_corner[j] < wave_min) { wave_min = wave_corner[j];}
      if(wave_corner[j] > wave_max) { wave_max = wave_corner[j];}
    }

    Area = find_area_quad(along_min, wave_min, along_corner, wave_corner);

    // estimate where the pixel overlaps in the cube
    // find the min and max values in the cube appropriate xcoord,ycoord or zcoord

    ia1 = (along_min - crval_along) / cdelt_along;
    ia2 = (along_max - crval_along) / cdelt_along;
    if (ia1 < 0){ ia1 = 0;}
    if (ia2 >= na){ ia2 = na -1;}

    MinW = (wave_min - crval3) / cdelt3;
    MaxW = (wave_max - crval3) / cdelt3;
    iz1= (int) MinW;
    iz2 = round(MaxW);

    if(iz1 < 0){ iz1 = 0;}
    if(iz2 < 0){ iz2 = iz1 +1;}
    if (iz2 >= nz){iz2 = nz - 1;}

     nplane = naxis1 * naxis2;
    // loop over possible overlapping cube pixels
    for(zz = iz1; zz < iz2+1; zz++){
      zcenter = zcoord[zz];
      istart = zz * nplane;
      for (aa = ia1; aa< ia2 + 1; aa++){
	cube_index = 0;
	if(instrument == 1) { // NIRSPec
	  cube_index = istart + aa * naxis1 + ss;  // ss = slice #
	} else {
	  cube_index = istart + ss * naxis1 + aa;   // yy = slice #
	}
	acenter = acoord[aa];
	area_overlap = sh_find_overlap(acenter, zcenter,
					      cdelt_along, cdelt3,
					      along_corner, wave_corner);


	if (area_overlap > 0.0) {
	  AreaRatio = area_overlap / Area;
	  fluxv[cube_index] = fluxv[cube_index] + (AreaRatio * pixel_flux[ipixel]);
	  weightv[cube_index] = weightv[cube_index] +	AreaRatio;
	  ifluxv[cube_index] = ifluxv[cube_index] + 1;
	  err = (AreaRatio * pixel_err[ipixel]) * (AreaRatio * pixel_err[ipixel]);
	  varv[cube_index] = varv[cube_index] + err;
	}
      }
    }

  }

  *spaxel_flux = fluxv;
  *spaxel_weight = weightv;
  *spaxel_var = varv;
  *spaxel_iflux = ifluxv;

  return 0;
}



// C extension SETUP

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


static PyObject *cube_wrapper_internal(PyObject *module, PyObject *args) {

  PyObject *result = NULL, *a1o, *a2o, *a3o, *a4o, *lam1o, *lam2o, *lam3o, *lam4o,
    *fluxo, *erro, *acoordo, *zcoordo;

  double crval_along, cdelt_along, crval3, cdelt3;
  int  nz, na, npt, naxis1, naxis2, ncube;
  int ss;
  int instrument;

  double *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  double *spaxel_iflux=NULL;

  int free_a1=0, free_a2=0, free_a3=0, free_a4 =0, free_lam1=0, free_lam2 =0, free_lam3=0, free_lam4=0;
  int free_acoord=0, free_zcoord=0, status=0;
  int free_flux=0, free_err=0;

  int n1, n2;
  PyArrayObject *a1=NULL, *a2=NULL, *a3=NULL, *a4=NULL, *lam1=NULL, *lam2=NULL;
  PyArrayObject *lam3=NULL, *lam4=NULL, *flux=NULL, *err=NULL, *acoord=NULL, *zcoord=NULL;

  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_var_arr=NULL;
  PyArrayObject *spaxel_iflux_arr=NULL;
  npy_intp npy_ncube = 0;

  if (!PyArg_ParseTuple(args, "iiiddddOOOOOOOOOOiOO:cube_wrapper_internal",
			&instrument, &naxis1, &naxis2, &crval_along, &cdelt_along,
			&crval3, &cdelt3,
			&a1o, &a2o, &a3o, &a4o,&lam1o, &lam2o, &lam3o, &lam4o,
			&acoordo, &zcoordo, &ss,  &fluxo, &erro)) {
    return NULL;
  }

  // check that input parameters are valid:
  if ((cdelt3 < 0) || (cdelt_along < 0)) {
    PyErr_SetString(PyExc_ValueError,
		    "'cdelt3' and 'cdelt_along' must be a strictly positive number.");
    return NULL;
  }

    // ensure we are working with numpy arrays and avoid creating new ones
    // if possible:
  if ((!(a1 = ensure_array(a1o, &free_a1))) ||
      (!(a2 = ensure_array(a2o, &free_a2))) ||
      (!(a3 = ensure_array(a3o, &free_a3))) ||
      (!(a4 = ensure_array(a4o, &free_a4))) ||
      (!(lam1 = ensure_array(lam1o, &free_lam1))) ||
      (!(lam2 = ensure_array(lam2o, &free_lam2))) ||
      (!(lam3 = ensure_array(lam3o, &free_lam3))) ||
      (!(lam4 = ensure_array(lam4o, &free_lam4))) ||
      (!(acoord = ensure_array(acoordo, &free_acoord))) ||
      (!(zcoord = ensure_array(zcoordo, &free_zcoord))) ||
      (!(flux = ensure_array(fluxo, &free_flux))) ||
      (!(err = ensure_array(erro, &free_err))) ) {
      goto cleanup;
    }

  npt = (int) PyArray_Size((PyObject *) flux);
  nz = (int) PyArray_Size((PyObject *) zcoord);
  na = (int) PyArray_Size((PyObject *) acoord);
  n1 = (int) PyArray_Size((PyObject *) a1);
  n2 = (int) PyArray_Size((PyObject *) lam1);

  if (n1 != npt || n2 != npt ) {
    PyErr_SetString(PyExc_ValueError,
		    "Input coordinate arrays of unequal size.");
    goto cleanup;
  }

  ncube = naxis1 * naxis2 * nz;

  if (ncube ==0) {
    // 0-length input arrays. Nothing to clip. Return 0-length arrays
    spaxel_flux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_DOUBLE, 0);
    if (!spaxel_flux_arr) goto fail;

    spaxel_weight_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_DOUBLE, 0);
    if (!spaxel_weight_arr) goto fail;

    spaxel_var_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_DOUBLE, 0);
    if (!spaxel_var_arr) goto fail;

    spaxel_iflux_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_DOUBLE, 0);
    if (!spaxel_iflux_arr) goto fail;

    result = Py_BuildValue("(NNNN)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
			   spaxel_iflux_arr);

    goto cleanup;
  }

  //______________________________________________________________________
  // Match the point cloud elements to the spaxels they fail within the roi
  //______________________________________________________________________
    status = match_detector_cube(instrument, naxis1, naxis2, nz, npt, ncube, na,
				 crval_along, cdelt_along, crval3, cdelt3,
				 (double *) PyArray_DATA(a1),
				 (double *) PyArray_DATA(a2),
				 (double *) PyArray_DATA(a3),
				 (double *) PyArray_DATA(a4),
				 (double *) PyArray_DATA(lam1),
				 (double *) PyArray_DATA(lam2),
				 (double *) PyArray_DATA(lam3),
				 (double *) PyArray_DATA(lam4),
				 (double *) PyArray_DATA(acoord),
				 (double *) PyArray_DATA(zcoord),
				 ss,
				 (double *) PyArray_DATA(flux),
				 (double *) PyArray_DATA(err),
				 &spaxel_flux, &spaxel_weight, &spaxel_var, &spaxel_iflux);

  if (status ) {
    goto fail;

  } else {
    // create return tuple:
    npy_ncube = (npy_intp) ncube;

    spaxel_flux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_DOUBLE, spaxel_flux);
    if (!spaxel_flux_arr) goto fail;
    spaxel_flux = NULL;

    spaxel_weight_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_DOUBLE, spaxel_weight);
    if (!spaxel_weight_arr) goto fail;
    spaxel_weight = NULL;

    spaxel_var_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_DOUBLE, spaxel_var);
    if (!spaxel_var_arr) goto fail;
    spaxel_var = NULL;

    spaxel_iflux_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_DOUBLE, spaxel_iflux);
    if (!spaxel_iflux_arr) goto fail;
    spaxel_iflux = NULL;

    PyArray_ENABLEFLAGS(spaxel_flux_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_weight_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_var_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_iflux_arr, NPY_ARRAY_OWNDATA);

    result = Py_BuildValue("(NNNN)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
			   spaxel_iflux_arr);

    goto cleanup;
  }

fail:
  Py_XDECREF(spaxel_flux_arr);
  Py_XDECREF(spaxel_weight_arr);
  Py_XDECREF(spaxel_var_arr);
  Py_XDECREF(spaxel_iflux_arr);

  free(spaxel_flux);
  free(spaxel_weight);
  free(spaxel_var);
  free(spaxel_iflux);

  if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_MemoryError,
		    "Unable to allocate memory for output arrays.");
  }

 cleanup:
  if (free_a1)Py_XDECREF(a1);
  if (free_a2)Py_XDECREF(a2);
  if (free_a3)Py_XDECREF(a3);
  if (free_a4)Py_XDECREF(a4);
  if (free_lam1)Py_XDECREF(lam1);
  if (free_lam2)Py_XDECREF(lam2);
  if (free_lam3)Py_XDECREF(lam3);
  if (free_lam4)Py_XDECREF(lam4);
  if (free_acoord) Py_XDECREF(acoord);
  if (free_zcoord) Py_XDECREF(zcoord);
  if (free_flux) Py_XDECREF(flux);
  if (free_err) Py_XDECREF(err);

  return result;
}

static PyMethodDef cube_methods[] =
{
    {
        "cube_wrapper_internal",
	cube_wrapper_internal,
	METH_VARARGS,
        "cube_wrapper_internal(put in doc string)"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cube_match_internal",             /* m_name */
    "find point cloud matches for each spaxel center",  /* m_doc */
    -1,                          /* m_size */
    cube_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_cube_match_internal(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
