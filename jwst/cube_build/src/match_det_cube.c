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
xcoord : numpy.ndarray
ycoord : numpy.ndarray
zcoord : numpy.ndarray
flux : numpy.ndarray
err : numpy.ndarray
coord1 : numpy.ndarray
coord2 : numpy.ndarray
wave : numpy.ndarray
rois_pixel : numpy.ndarray
roiw_pixel : numpy.ndarray
scalerad_pixel : numpy.ndarray
zcdelt3: numpy.ndarray

nplane : int
nwave : int
ncube : int (this is a huge number does it need to int64?)
npt : int  

cdelt1 : double
cdelt2 : double


Returns
-------

spaxel_flux : numpy.ndarray
spaxel_weight : numpy.ndarray
spaxel_iflux : numpy.ndarray
spaxel_var : numpy.ndarray

*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_match_det_cube_numpy_api    //WHAT IS THIS AND WHERE IS IT USED???
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>



int mem_alloc(int nelem, double **fluxv, double **weightv, double **varv,  double **ifluxv) {

    double *f, *w, *v, *i;
    const char *msg = "Couldn't allocate memory for output arrays.";

    // flux:
    f = (double*)malloc(nelem* sizeof(double));
    if (f) {
        *fluxv = f;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // weight:
    w = (double*)malloc(nelem* sizeof(double));
    if (w) {
        *weightv = w;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // variance:
    v = (double*)malloc(nelem* sizeof(double));
    if (v) {
        *varv = v;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }
    // iflux
    i = (double*)malloc(nelem* sizeof(double));
    if (i) {
        *ifluxv = i;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    return 0;
}



/*
return values: spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux
*/
int match_point_emsm(double *xc, double *yc, double *zc,
		     double *coord1, double *coord2, double *wave,
		     double *flux, double *err,
		     double *rois_pixel, double *roiw_pixel, double *scalerad_pixel,
		     double *zcdelt3,
		     int nx, int ny, int nwave, int ncube, int npt,
		     double cdelt1, double cdelt2,
		     double **spaxel_flux, double **spaxel_weight, double **spaxel_var,
		     double **spaxel_iflux) {

    double *fluxv = NULL, *weightv=NULL, *varv=NULL ;  // vectors for spaxel 
    double *ifluxv = NULL;  // vector for spaxel

    // allocate memory to hold output 
    //if (mem_alloc(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;
    int status = (mem_alloc(ncube, &fluxv, &weightv, &varv, &ifluxv));
    printf("Return array from allocating memory %i \n", status);
    
    // Set all data to zero
    for (int i = 0; i < ncube; i++){
      varv[i] = 0.0;
      fluxv[i] = 0.0;
      ifluxv[i] = 0.0;
      weightv[i] = 0.0;
    }
    
    // loop over each point cloud member and find which roi spaxels it is found

    for (int k = 0; k < npt; k++) {
       // Search wave and find match
      int iwstart = -1;
      int iwend = -1;
      int ii = 0;
      int done_search_w = 0;
      float check_lower = zc[0] - roiw_pixel[k];
      float check_upper = zc[nwave-1] + roiw_pixel[k];
      if(wave[k] < check_lower || wave[k] > check_upper){
	printf(" point cloud member wavelength range outside IFU cube %f %f %f \n",check_lower, check_upper, wave[k]);
	  }
      while (ii < nwave && done_search_w == 0) {
	float wdiff = fabs(zc[ii] - wave[k]);
	if(wdiff <= roiw_pixel[k]){
	  if (iwstart == -1){
	    iwstart = ii;
	  }
	} else{
	  if(iwstart != -1 && iwend ==-1){
	    iwend = ii;
	    done_search_w = 1;
	  }
	}
	ii = ii + 1;
      }
      // catch the case of iwstart near nwave and becomes = nwave before iwend can be set.
      if(iwstart !=-1 && iwend == -1){
      	iwend = nwave;
      	done_search_w = 1;
      }

      
      // Search xcenters and find match
      int ixstart = -1;
      int ixend = -1;
      ii = 0;
      int done_search_x = 0;
      while (ii < nx && done_search_x == 0) {
	double xdiff = fabs(xc[ii] - coord1[k]);
	if(xdiff <= rois_pixel[k]){
	  if (ixstart == -1){
	    ixstart = ii;
	  }
	} else{
	  if(ixstart != -1 && ixend ==-1){
	    ixend = ii;
	    done_search_x = 1;
	  }
	}
	ii = ii + 1;
      }
      // catch the case of ixstart near nx and becomes = nx before ixend can be set.
      if(ixstart !=-1 && ixend == -1){
	ixend = nx;
	done_search_x = 1;
      }
      
       // Search ycenters and find match
      int iystart = -1;
      int iyend = -1;
      ii = 0;
      int done_search_y = 0;
      while (ii < ny && done_search_y == 0) {
	double ydiff = fabs(yc[ii] - coord2[k]);
	if(ydiff <= rois_pixel[k]){
	  if (iystart == -1){
	    iystart = ii;
	  }
	} else{
	  if(iystart != -1 && iyend ==-1){
	    iyend = ii;
	    done_search_y = 1;
	  }
	}
	ii = ii + 1;
      }
      // catch the case of iystart near ny and becomes = ny before iyend can be set.
      if(iystart !=-1 && iyend == -1){
	iyend = ny;
	done_search_y = 1;
      }

      // set up the values for fluxv, weightv, ifluxv, varv
      int nxy = nx * ny;
      if(done_search_x == 1 && done_search_y ==1 && done_search_w ==1){
	// The search above for x,y  was a crude search - now narrow the search using the distance between
	// the spaxel center and point cloud
	for (int ix = ixstart; ix< ixend; ix ++){
	  for ( int iy = iystart; iy < iyend; iy ++){
	    double ydist = fabs(yc[iy] - coord2[k]);
	    double xdist = fabs(xc[ix] - coord1[k]);
	    double radius = sqrt( xdist*xdist + ydist*ydist);

	    if (radius <= rois_pixel[k]){
	      // Find the index for this in spatial plane
 	      int index_xy = iy* nx + ix;
	      for (int iw = iwstart; iw< iwend; iw++){
		int index_cube = iw*nxy + index_xy;
		if(index_cube > ncube){
		  printf(" Index Cube > ncube \n");
		}
		double d1 = xdist/cdelt1;
		double d2 = ydist/cdelt2;
		double dxy = (d1 * d1) + (d2 * d2);
		double d3 = (wave[k] - zc[iw])/ zcdelt3[iw];
		double wdist = fabs( wave[k] - zc[iw]);
		double d32 = d3 * d3;
		double w = d32  +  dxy;
		double wn = -w/(scalerad_pixel[k]/cdelt1);
		double ww = exp(wn);

		double weighted_flux =  flux[k]* ww;
		double weighted_var = (err[k]* ww) * (err[k]*ww);
		fluxv[index_cube] = fluxv[index_cube] + weighted_flux;
		weightv[index_cube] = weightv[index_cube] + ww;
		varv[index_cube] = varv[index_cube] + weighted_var;
		ifluxv[index_cube] = ifluxv[index_cube] +1.0;

		if( ix == 24 && iy == 24 && iw ==0) {
		  printf("found element %i %i %f %f %f  %f %f \n", k, index_cube, err[k], flux[k], ifluxv[index_cube], varv[index_cube], fluxv[index_cube]);
		}
	      }
	    }
	  } // end loop over iy
	} // end loop over ix

      } // end done_search_x, done_search_y, done_search_w
    } // end loop over point cloud
    

    // assign output values:

    *spaxel_flux = fluxv;
    *spaxel_weight = weightv;
    *spaxel_var = varv;
    *spaxel_iflux = ifluxv;

    free(fluxv);
    free(weightv);
    free(varv);
    free(ifluxv);
    
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


static PyObject *point_emsm(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xco, *yco, *zco, *fluxo, *erro, *coord1o, *coord2o, *waveo;
  PyObject *rois_pixelo, *roiw_pixelo, *scalerad_pixelo, *zcdelt3o;
  
  double cdelt1, cdelt2;
  int  nwave, npt, nxx, nyy, ncube;
  double *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  double *spaxel_iflux=NULL;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_rois_pixel=0, free_roiw_pixel=0, free_scalerad_pixel=0, free_flux=0, free_err=0, free_zcdelt3=0;
  
  PyArrayObject *xc, *yc, *zc, *flux, *err, *coord1, *coord2, *wave, *rois_pixel, *roiw_pixel, *scalerad_pixel;
  PyArrayObject *zcdelt3;
  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_var_arr=NULL;
  PyArrayObject *spaxel_iflux_arr=NULL; 
  npy_intp npy_ncube = 0;

  int  ny,nz;
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOddiiiii:point_emsm",
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro,
			&rois_pixelo, &roiw_pixelo, &scalerad_pixelo,&zcdelt3o,
			&cdelt1, &cdelt2, &nxx, &nyy, &nwave, &npt, &ncube)) {
    return NULL;
  }


  // check that input parameters are valid:

  if ((cdelt1 < 0) || (cdelt2 < 0)) {
    PyErr_SetString(PyExc_ValueError,
		    "'cdelt1' and 'cdelt2' must be a strictly positive number.");
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
      (!(err = ensure_array(erro, &free_err))) ||
      (!(rois_pixel = ensure_array(rois_pixelo, &free_rois_pixel))) ||
      (!(roiw_pixel = ensure_array(roiw_pixelo, &free_roiw_pixel))) ||
      (!(zcdelt3 = ensure_array(zcdelt3o, &free_zcdelt3))) ||
      (!(scalerad_pixel = ensure_array(scalerad_pixelo, &free_scalerad_pixel))) )
    {
      goto cleanup;
    }

  npt = (int) PyArray_Size((PyObject *) coord1);
  ny = (int) PyArray_Size((PyObject *) coord2);
  nz = (int) PyArray_Size((PyObject *) wave);
  if (npt != ny) {
    PyErr_SetString(PyExc_ValueError,
		    "Input coordinate arrays of unequal size.");
    goto cleanup;
  }

  //nxx = (int) PyArray_Size((PyObject *) xc);
  //nyy = (int) PyArray_Size((PyObject *) yc);
  //nwave = (int) PyArray_Size((PyObject *) zc);

  //ncube = nxx * nyy * nwave;
  
  printf(" sizes %i %i %i %i %i \n ", nxx, nyy, nwave, npt, ncube);

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

    result = Py_BuildValue("(OOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr, spaxel_iflux_arr);
	
    goto cleanup;
  }

  status = match_point_emsm((double *) PyArray_DATA(xc),
			    (double *) PyArray_DATA(yc),
			    (double *) PyArray_DATA(zc),
			    (double *) PyArray_DATA(coord1),
			    (double *) PyArray_DATA(coord2),
			    (double *) PyArray_DATA(wave),			      
			    (double *) PyArray_DATA(flux),
			    (double *) PyArray_DATA(err),
			    (double *) PyArray_DATA(rois_pixel),
			    (double *) PyArray_DATA(roiw_pixel),
			    (double *) PyArray_DATA(scalerad_pixel),
			    (double *) PyArray_DATA(zcdelt3),
			    nxx, nyy, nwave, ncube, npt, cdelt1, cdelt2,
			    &spaxel_flux, &spaxel_weight, &spaxel_var, &spaxel_iflux);

  if (status) {
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
    
    result = Py_BuildValue("(OOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr, spaxel_iflux_arr);
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
  if (free_xc) Py_XDECREF(xc);
  if (free_yc) Py_XDECREF(yc);
  if (free_zc) Py_XDECREF(zc);
  if (free_coord1) Py_XDECREF(coord1);
  if (free_coord2) Py_XDECREF(coord2);
  if (free_wave) Py_XDECREF(wave);
  if (free_flux) Py_XDECREF(flux);
  if (free_err) Py_XDECREF(err);
  if (free_rois_pixel) Py_XDECREF(rois_pixel);
  if (free_roiw_pixel) Py_XDECREF(roiw_pixel);
  if (free_scalerad_pixel) Py_XDECREF(scalerad_pixel);
  if (free_zcdelt3) Py_XDECREF(zcdelt3);
  return result;
}


static PyMethodDef match_det_cube_methods[] =
{
    {
        "point_emsm",
	point_emsm,
	METH_VARARGS,
        "point_emsm(put in doc string)"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
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
