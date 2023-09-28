/*
The detector pixels are represented by a 'point cloud' on the sky. The IFU cube is
represented by a 3-D regular grid. This module finds the point cloud members contained
in a region centered on the center of the cube spaxel. The size of the spaxel is spatial
coordinates is cdetl1 and cdelt2, while the wavelength size is zcdelt3.
This module uses the modified shephard weighting method (emsm if weight_type =0 or msm if weight_type =1)
to determine how to  weight each point cloud member in the spaxel.

Main function for Python: cube_wrapper

Python signature: result = cube_wrapper(instrument, flag_dq_plane, weight_type, start_region, end_region,
                                        overlap_partial, overlap_full,
                                        xcoord, ycoord, zcoord,
                                        coord1, coord2, wave, flux, err, slice_no,
                                        rois_pixel, roiw_pixel, scalerad_pixel
					weight_pixel, softrad_pixel,cdelt3_normal,
                                        roiw_ave, cdelt1, cdelt2)
provide more details

The output of this function is a tuple of 5 arrays:(spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq)
example output

Parameters
----------
instrument : int
    0 = MIRI, 1 = NIRSPEC. Used for set the dq plane
flag_dq_plane : int
   0 do set the DQ plane based on FOV, but set all values =0
   1 set the DQ plane based on the FOV
weight_type : int
   0: use emsm weighting
   1: use msm weighting
start_region : int
    starting slice number for detector region used in dq flagging
end_region: int
    ending slice number for detector region used in dq flagging
overlap_partial : int
    a dq flag indicating that only a portion of the spaxel is overlapped by a mapped detector pixel
overlap_full : int
    a dq flag indicating that the entire spaxel is overlapped by the mapped detector pixel
xcoord : double array
   size of naxis1. This array holds the center x axis values of the ifu cube
ycoord : double array
   size of naxis2. This array holds the center y axis values of the ifu cube
zcoord : double array
   size of naxis3. This array holds the center x axis values of the ifu cube
flux : double array
   size: point cloud elements. Flux of each point cloud member
err : double array
   size: point cloud elements. err of each point cloud member
slice_no: int
   slice number of point cloud member to be in dq flagging
coord1 : double array
   size: point cloud elements. Naxis 1 coordinate of point cloud member (xi)
coord2 : double array
   size: point cloud elements. Naxis 2 coordinate of point cloud member (eta)
wave : double array
   size: point cloud elements. Wavelength of each point cloud member
rois_pixel : double array
   size: point cloud elements. Roi in spatial dimension to use for point cloud member
roiw_pixel : double array
   size: point cloud elements. Roi in wavelength dimension to use point cloud member
scalerad_pixel : double array
   size: point cloud elements. MSM weight parameter to use for point cloud member
zcdelt3: double array
   size: point cloud elements. Spectral scale to use at wavelength of point cloud member
roiw_ave : double
   Average roiw for all the wavelength planes. Used in dq flagging
cdelt1 : double
   Naxis 1 scale for cube
cdelt2 : double
   Naxis 2 scale for cube


Returns
-------
spaxel_flux : numpy.ndarray
  IFU spaxel cflux
spaxel_weight : numpy.ndarray
  IFU spaxel weight
spaxel_iflux : numpy.ndarray
  IFU spaxel weight map (number of overlaps)
spaxel_var : numpy.ndarray
  IFU spaxel error
spaxel_dq : numpy.ndarray
  IFU spaxel dq
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_match_sky_pointcloud_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// routines used from cube_utils.c

extern int alloc_flux_arrays(int nelem, double **fluxv, double **weightv, double **varv,  double **ifluxv);


extern int dq_miri(int start_region, int end_region, int overlap_partial, int overlap_full,
		   int nx, int ny, int nz,
		   double cdelt1, double cdelt2, double roiw_ave,
		   double *xc, double *yc, double *zc,
		   double *coord1, double *coord2, double *wave,
		   double *sliceno,
		   long ncube, long npt,
		   int **spaxel_dq);

extern int dq_nirspec(int overlap_partial,
		      int nx, int ny, int nz,
		      double cdelt1, double cdelt2, double roiw_ave,
		      double *xc, double *yc, double *zc,
		      double *coord1, double *coord2, double *wave,
		      double *sliceno,
		      long ncube, long npt,
		      int **spaxel_dq);

extern int set_dqplane_to_zero(int ncube, int **spaxel_dq);

// Match point cloud to sky and determine the weighting to assign to each point cloud  member
// to matched spaxel based on ROI - weighting type - emsm

// return values: spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux

int match_point_emsm(double *xc, double *yc, double *zc,
		     double *coord1, double *coord2, double *wave,
		     double *flux, double *err,
		     double *rois_pixel, double *roiw_pixel, double *scalerad_pixel,
		     double *zcdelt3,
		     int nx, int ny, int nwave, int ncube, int npt,
		     double cdelt1, double cdelt2,
		     double **spaxel_flux, double **spaxel_weight, double **spaxel_var,
		     double **spaxel_iflux) {


  double *fluxv=NULL, *weightv=NULL, *varv=NULL, *ifluxv=NULL;  // vector for spaxel

  int k, iwstart, iwend, ixstart,  ixend, iystart, iyend;
  int ii, nxy, ix, iy, iw, index_xy, index_cube;
  int done_search_w, done_search_y, done_search_x;
  double wdiff, ydiff, xdiff, ydist, xdist, radius;
  double d1, d2, dxy, d3, d32, w, wn, ww, weighted_flux, weighted_var;

  // allocate memory to hold output
  if (alloc_flux_arrays(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;


    // loop over each point cloud member and find which roi spaxels it is found

  for (k = 0; k < npt; k++) {
    // Search wave and find match
    iwstart = -1;
    iwend = -1;
    ii = 0;
    done_search_w = 0;

    while (ii < nwave && done_search_w == 0) {
      wdiff = fabs(zc[ii] - wave[k]);
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
    ixstart = -1;
    ixend = -1;
    ii = 0;
    done_search_x = 0;

    while (ii < nx && done_search_x == 0) {
      xdiff = fabs(xc[ii] - coord1[k]);
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
    iystart = -1;
    iyend = -1;
    ii = 0;
    done_search_y = 0;
    while (ii < ny && done_search_y == 0) {
      ydiff = fabs(yc[ii] - coord2[k]);
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
    nxy = nx * ny;
    if(done_search_x == 1 && done_search_y ==1 && done_search_w ==1){
      // The search above for x,y  was a crude search - now narrow the search using the distance between
      // the spaxel center and point cloud

      for (ix = ixstart; ix< ixend; ix ++){

	for ( iy = iystart; iy < iyend; iy ++){
	  ydist = fabs(yc[iy] - coord2[k]);
	  xdist = fabs(xc[ix] - coord1[k]);
	  radius = sqrt( xdist*xdist + ydist*ydist);

	  if (radius <= rois_pixel[k]){
	    // Find the index for this in spatial plane
	    index_xy = iy* nx + ix;
	    for (iw = iwstart; iw< iwend; iw++){
	      index_cube = iw*nxy + index_xy;

	      d1 = xdist/cdelt1;
	      d2 = ydist/cdelt2;
	      dxy = (d1 * d1) + (d2 * d2);
	      d3 = (wave[k] - zc[iw])/ zcdelt3[iw];
	      d32 = d3 * d3;
	      w = d32  +  dxy;
	      wn = -w/(scalerad_pixel[k]/cdelt1);
	      ww = exp(wn);

	      weighted_flux =  flux[k]* ww;
	      weighted_var = (err[k]* ww) * (err[k]*ww);
	      fluxv[index_cube] = fluxv[index_cube] + weighted_flux;
	      weightv[index_cube] = weightv[index_cube] + ww;
	      varv[index_cube] = varv[index_cube] + weighted_var;
	      ifluxv[index_cube] = ifluxv[index_cube] +1.0;
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

  return 0;
}



// Match point cloud to sky and determine the weighting to assign to each point cloud  member
// to matched spaxel based on ROI - weighting type - msm
//return values: spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux

int match_point_msm(double *xc, double *yc, double *zc,
		    double *coord1, double *coord2, double *wave,
		    double *flux, double *err,
		    double *rois_pixel, double *roiw_pixel,
		    double *weight_pixel, double *softrad_pixel,
		    double *zcdelt3,
		    int nx, int ny, int nwave, int ncube, int npt,
		    double cdelt1, double cdelt2,
		    double **spaxel_flux, double **spaxel_weight, double **spaxel_var,
		    double **spaxel_iflux) {


  double *fluxv=NULL, *weightv=NULL, *varv=NULL, *ifluxv=NULL;  // vector for spaxel

  int k;
  int iwstart, iwend, ixstart, ixend, iystart, iyend;
  int ii, nxy, iw, ix, iy, index_xy, index_cube;
  int done_search_w, done_search_x, done_search_y;
  double wdiff, xdiff, ydiff, radius, ydist, xdist;
  double d1, d2, dxy, d3, d32, w, wn, ww;
  double weighted_flux, weighted_var;

  // allocate memory to hold output
  if (alloc_flux_arrays(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;

  // loop over each point cloud member and find which roi spaxels it is found

  for ( k = 0; k < npt; k++) {
    // Search wave and find match
    iwstart = -1;
    iwend = -1;
    ii = 0;
    done_search_w = 0;

    while (ii < nwave && done_search_w == 0) {
	wdiff = fabs(zc[ii] - wave[k]);
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
      ixstart = -1;
      ixend = -1;
      ii = 0;
      done_search_x = 0;
      while (ii < nx && done_search_x == 0) {
	xdiff = fabs(xc[ii] - coord1[k]);
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
      iystart = -1;
      iyend = -1;
      ii = 0;
      done_search_y = 0;
      while (ii < ny && done_search_y == 0) {
	ydiff = fabs(yc[ii] - coord2[k]);
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
      nxy = nx * ny;
      if(done_search_x == 1 && done_search_y ==1 && done_search_w ==1){
	// The search above for x,y  was a crude search - now narrow the search using the distance between
	// the spaxel center and point cloud
	for (ix = ixstart; ix< ixend; ix ++){
	  for ( iy = iystart; iy < iyend; iy ++){
	    ydist = fabs(yc[iy] - coord2[k]);
	    xdist = fabs(xc[ix] - coord1[k]);
	    radius = sqrt( xdist*xdist + ydist*ydist);

	    if (radius <= rois_pixel[k]){
	      // Find the index for this in spatial plane
 	      index_xy = iy* nx + ix;
	      for (iw = iwstart; iw< iwend; iw++){
		index_cube = iw*nxy + index_xy;

		d1 = xdist/cdelt1;
		d2 = ydist/cdelt2;
		dxy = (d1 * d1) + (d2 * d2);
		d3 = (wave[k] - zc[iw])/ zcdelt3[iw];

		d32 = d3 * d3;
		w = d32  +  dxy;
		wn = pow(sqrt(w), weight_pixel[k]);
		if( wn < softrad_pixel[k]){
		  wn = softrad_pixel[k];
		}

		ww = 1.0/wn;
		weighted_flux =  flux[k]* ww;
		weighted_var = (err[k]* ww) * (err[k]*ww);
		fluxv[index_cube] = fluxv[index_cube] + weighted_flux;
		weightv[index_cube] = weightv[index_cube] + ww;
		varv[index_cube] = varv[index_cube] + weighted_var;
		ifluxv[index_cube] = ifluxv[index_cube] +1.0;

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


static PyObject *cube_wrapper(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xco, *yco, *zco, *fluxo, *erro, *coord1o, *coord2o, *waveo, *slicenoo;
  PyObject *rois_pixelo, *roiw_pixelo, *scalerad_pixelo, *zcdelt3o, *softrad_pixelo, *weight_pixelo;

  double cdelt1, cdelt2, roiw_ave;
  int  nwave, npt, nxx, nyy, ncube;

  int instrument, flag_dq_plane,start_region, end_region, overlap_partial, overlap_full, weight_type;
  double *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  double *spaxel_iflux=NULL;
  int *spaxel_dq=NULL;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_rois_pixel=0, free_roiw_pixel=0, free_scalerad_pixel=0, free_flux=0, free_err=0, free_zcdelt3=0;
  int free_sliceno=0, free_softrad_pixel, free_weight_pixel=0;

  PyArrayObject *xc, *yc, *zc, *flux, *err, *coord1, *coord2, *wave, *rois_pixel, *roiw_pixel, *scalerad_pixel;
  PyArrayObject *zcdelt3, *sliceno, *softrad_pixel, *weight_pixel;
  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_var_arr=NULL;
  PyArrayObject *spaxel_iflux_arr=NULL, *spaxel_dq_arr=NULL;
  npy_intp npy_ncube = 0;

  int  ny,nz;

  if (!PyArg_ParseTuple(args, "iiiiiiiOOOOOOOOOOOOOOOddd:cube_wrapper",
			&instrument, &flag_dq_plane, &weight_type,  &start_region, &end_region, &overlap_partial, &overlap_full,
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro, &slicenoo,
			&rois_pixelo, &roiw_pixelo, &scalerad_pixelo, &weight_pixelo, &softrad_pixelo, &zcdelt3o, &roiw_ave,
			&cdelt1, &cdelt2)) {
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
      (!(sliceno = ensure_array(slicenoo, &free_sliceno))) ||
      (!(rois_pixel = ensure_array(rois_pixelo, &free_rois_pixel))) ||
      (!(roiw_pixel = ensure_array(roiw_pixelo, &free_roiw_pixel))) ||
      (!(zcdelt3 = ensure_array(zcdelt3o, &free_zcdelt3))) ||
      (!(scalerad_pixel = ensure_array(scalerad_pixelo, &free_scalerad_pixel))) ||
      (!(softrad_pixel = ensure_array(softrad_pixelo, &free_softrad_pixel))) ||
      (!(weight_pixel = ensure_array(weight_pixelo, &free_weight_pixel))) )
    {
      goto cleanup;

    }

  npt = (int) PyArray_Size((PyObject *) coord1);
  ny = (int) PyArray_Size((PyObject *) coord2);
  nz = (int) PyArray_Size((PyObject *) wave);
  if (ny != npt || nz != npt ) {
    PyErr_SetString(PyExc_ValueError,
		    "Input coordinate arrays of unequal size.");
    goto cleanup;

  }

  nxx = (int) PyArray_Size((PyObject *) xc);
  nyy = (int) PyArray_Size((PyObject *) yc);
  nwave = (int) PyArray_Size((PyObject *) zc);

  ncube = nxx * nyy * nwave;

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

    spaxel_dq_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_INT, 0);
    if (!spaxel_dq_arr) goto fail;

    result = Py_BuildValue("(NNNNN)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
			   spaxel_iflux_arr, spaxel_dq_arr);

    goto cleanup;

  }

  //______________________________________________________________________
  // if flag_dq_plane = 1, Set up the dq plane
  //______________________________________________________________________
  int status1 = 0;
  if(flag_dq_plane){
    if (instrument == 0){
      status1 = dq_miri(start_region, end_region,overlap_partial, overlap_full,
			nxx,nyy,nwave,
			cdelt1, cdelt2, roiw_ave,
			(double *) PyArray_DATA(xc),
			(double *) PyArray_DATA(yc),
			(double *) PyArray_DATA(zc),
			(double *) PyArray_DATA(coord1),
			(double *) PyArray_DATA(coord2),
			(double *) PyArray_DATA(wave),
			(double *) PyArray_DATA(sliceno),
			ncube, npt,
			&spaxel_dq);

    } else{
      status1 = dq_nirspec(overlap_partial,
			   nxx,nyy,nwave,
			   cdelt1, cdelt2, roiw_ave,
			   (double *) PyArray_DATA(xc),
			   (double *) PyArray_DATA(yc),
			   (double *) PyArray_DATA(zc),
			   (double *) PyArray_DATA(coord1),
			   (double *) PyArray_DATA(coord2),
			   (double *) PyArray_DATA(wave),
			   (double *) PyArray_DATA(sliceno),
			   ncube, npt,
			   &spaxel_dq);
    }
  } else{ // set dq plane to 0
    status1 = set_dqplane_to_zero(ncube, &spaxel_dq);

  }

  //______________________________________________________________________
  // Match the point cloud elements to the spaxels they fail within the roi
  //______________________________________________________________________
  if(weight_type ==0){
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
  } else{
    status = match_point_msm((double *) PyArray_DATA(xc),
			     (double *) PyArray_DATA(yc),
			     (double *) PyArray_DATA(zc),
			     (double *) PyArray_DATA(coord1),
			     (double *) PyArray_DATA(coord2),
			     (double *) PyArray_DATA(wave),
			     (double *) PyArray_DATA(flux),
			     (double *) PyArray_DATA(err),
			     (double *) PyArray_DATA(rois_pixel),
			     (double *) PyArray_DATA(roiw_pixel),
			     (double *) PyArray_DATA(weight_pixel),
			     (double *) PyArray_DATA(softrad_pixel),
			     (double *) PyArray_DATA(zcdelt3),
			     nxx, nyy, nwave, ncube, npt, cdelt1, cdelt2,
			     &spaxel_flux, &spaxel_weight, &spaxel_var, &spaxel_iflux);
  }


  if (status || status1) {
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

    spaxel_dq_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_INT, spaxel_dq);
    if (!spaxel_dq_arr) goto fail;
    spaxel_dq = NULL;

    PyArray_ENABLEFLAGS(spaxel_flux_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_weight_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_var_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_iflux_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS(spaxel_dq_arr, NPY_ARRAY_OWNDATA);
    result = Py_BuildValue("(NNNNN)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
			   spaxel_iflux_arr, spaxel_dq_arr);

    goto cleanup;

  }

 fail:
  Py_XDECREF(spaxel_flux_arr);
  Py_XDECREF(spaxel_weight_arr);
  Py_XDECREF(spaxel_var_arr);
  Py_XDECREF(spaxel_iflux_arr);
  Py_XDECREF(spaxel_dq_arr);

  free(spaxel_flux);
  free(spaxel_weight);
  free(spaxel_var);
  free(spaxel_iflux);
  free(spaxel_dq);

  if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_MemoryError,
		    "Unable to allocate memory for output arrays.");
  }

 cleanup:
  if (free_xc)Py_XDECREF(xc);
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
  if (free_softrad_pixel) Py_XDECREF(softrad_pixel);
  if (free_weight_pixel) Py_XDECREF(weight_pixel);
  if (free_zcdelt3) Py_XDECREF(zcdelt3);
  if (free_sliceno) Py_XDECREF(sliceno);

  return result;
}


static PyMethodDef cube_methods[] =
{
    {
        "cube_wrapper",
	cube_wrapper,
	METH_VARARGS,
        "cube_wrapper(put in doc string)"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cube_match_sky_pointcloud",  /* m_name */
    "find point cloud matches for each spaxel center",  /* m_doc */
    -1,                          /* m_size */
    cube_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_cube_match_sky_pointcloud(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
