/*
The detector pixels are represented by a 'point cloud' on the sky. The IFU cube is
represented by a 3-D regular grid. This module finds the point cloud members contained
in a region centered on the center of the cube spaxel. The size of the spaxel is spatial
coordinates is cdetl1 and cdelt2, while the wavelength size is zcdelt3.
This module uses the modified shephard weighting method (emsm if weight_type =0 or msm if weight_type =1)
to determine how to  weight each point cloud member in the spaxel.

Main function for Python: cube_wrapper

Python signature:result = cube_wrapper(instrument, flag_dq_plane, weight_type, start_region, end_region,
                                       self.overlap_partial, self.overlap_full,
                                       self.xcoord, self.ycoord, self.zcoord,
                                       coord1, coord2, wave, flux, err, slice_no,
                                       rois_pixel, roiw_pixel, scalerad_pixel,
                                       weight_pixel, softrad_pixel,
                                       self.cdelt3_normal,
                                       spaxel_flux, spaxel_weight, spaxel_var,
                                       spaxel_iflux, spaxel_dq,
                                       roiw_ave, self.cdelt1, self.cdelt2)

The arrays, spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq, are passed in a zero filled arrays
and filled in by this routine. 

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


extern int dq_miri(int start_region, int end_region, int overlap_partial, int overlap_full,
		   int nx, int ny, int nz,
		   double cdelt1, double cdelt2, double roiw_ave,
		   double *xc, double *yc, double *zc,
		   double *coord1, double *coord2, double *wave,
		   double *sliceno,
		   long ncube, int npt,
		   int *spaxel_dq);

extern int dq_nirspec(int overlap_partial,
		      int nx, int ny, int nz,
		      double cdelt1, double cdelt2, double roiw_ave,
		      double *xc, double *yc, double *zc,
		      double *coord1, double *coord2, double *wave,
		      double *sliceno,
		      long ncube, int npt,
		      int *spaxel_dq);


// Match point cloud to sky and determine the weighting to assign to each point cloud  member
// to matched spaxel based on ROI - weighting type - emsm

int match_point_emsm(double *xc, double *yc, double *zc,
		     double *coord1, double *coord2, double *wave,
		     double *flux, double *err,
		     double *rois_pixel, double *roiw_pixel, double *scalerad_pixel,
		     double *zcdelt3,
		     int nx, int ny, int nwave, int ncube, int npt,
		     double cdelt1, double cdelt2,
		     double *spaxel_flux, double *spaxel_weight, double *spaxel_var,
		     double *spaxel_iflux) {


  int k, iwstart, iwend, ixstart,  ixend, iystart, iyend;
  int ii, nxy, ix, iy, iw, index_xy, index_cube;
  int done_search_w, done_search_y, done_search_x;
  double wdiff, ydiff, xdiff, ydist, xdist, radius;
  double d1, d2, dxy, d3, d32, w, wn, ww, weighted_flux, weighted_var;

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
	      spaxel_flux[index_cube] = spaxel_flux[index_cube] + weighted_flux;
	      spaxel_weight[index_cube] = spaxel_weight[index_cube] + ww;
	      spaxel_var[index_cube] = spaxel_var[index_cube] + weighted_var;
	      spaxel_iflux[index_cube] = spaxel_iflux[index_cube] +1.0;
	    }
	  }
	} // end loop over iy
      } // end loop over ix
    } // end done_search_x, done_search_y, done_search_w
  } // end loop over point cloud

  // assign output values:

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
		    double *spaxel_flux, double *spaxel_weight, double *spaxel_var,
		    double *spaxel_iflux) {


  int k;
  int iwstart, iwend, ixstart, ixend, iystart, iyend;
  int ii, nxy, iw, ix, iy, index_xy, index_cube;
  int done_search_w, done_search_x, done_search_y;
  double wdiff, xdiff, ydiff, radius, ydist, xdist;
  double d1, d2, dxy, d3, d32, w, wn, ww;
  double weighted_flux, weighted_var;


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
		spaxel_flux[index_cube] = spaxel_flux[index_cube] + weighted_flux;
		spaxel_weight[index_cube] = spaxel_weight[index_cube] + ww;
		spaxel_var[index_cube] = spaxel_var[index_cube] + weighted_var;
		spaxel_iflux[index_cube] = spaxel_iflux[index_cube] +1.0;
	      }
	    }
	  } // end loop over iy
	} // end loop over ix

      } // end done_search_x, done_search_y, done_search_w
    } // end loop over point cloud

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

PyArrayObject * ensure_array_int(PyObject *obj, int *is_copy) {
  if (PyArray_CheckExact(obj) &&
      PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
      PyArray_TYPE((PyArrayObject *) obj) == NPY_UINT32) {
    *is_copy = 0;
    return (PyArrayObject *) obj;
  } else {
    *is_copy = 1;
    return (PyArrayObject *) PyArray_FromAny(obj, PyArray_DescrFromType(NPY_UINT32), 0, 0,
					     NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST, NULL
					     );
  }
}

static PyObject *cube_wrapper(PyObject *module, PyObject *args) {
  
  PyObject *xco, *yco, *zco, *fluxo, *erro, *coord1o, *coord2o, *waveo, *slicenoo;
  PyObject *rois_pixelo, *roiw_pixelo, *scalerad_pixelo, *zcdelt3o, *softrad_pixelo, *weight_pixelo;
  PyObject *spaxel_fluxo, *spaxel_weighto, *spaxel_varo, *spaxel_ifluxo;
  PyObject *spaxel_dqo;
  
  double cdelt1, cdelt2, roiw_ave;
  int  nwave, npt, nxx, nyy, ncube;

  int instrument, flag_dq_plane,start_region, end_region, overlap_partial, overlap_full, weight_type;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_rois_pixel=0, free_roiw_pixel=0, free_scalerad_pixel=0, free_flux=0, free_err=0, free_zcdelt3=0;
  int free_sliceno=0, free_softrad_pixel, free_weight_pixel=0;
  int free_spaxel_flux=0, free_spaxel_weight=0, free_spaxel_var=0, free_spaxel_iflux=0, free_spaxel_dq=0;

  PyArrayObject *xc, *yc, *zc, *flux, *err, *coord1, *coord2, *wave, *rois_pixel, *roiw_pixel, *scalerad_pixel;
  PyArrayObject *zcdelt3, *sliceno, *softrad_pixel, *weight_pixel;
  PyArrayObject *spaxel_flux, *spaxel_weight, *spaxel_var, *spaxel_iflux, *spaxel_dq;
  
  const int max_size_error = 80;
  char error[max_size_error] = "None";
  int flag_error  = 0;

  int  ny,nz;

  if (!PyArg_ParseTuple(args, "iiiiiiiOOOOOOOOOOOOOOOOOOOOddd:cube_wrapper",
			&instrument, &flag_dq_plane, &weight_type,  &start_region, &end_region, &overlap_partial, &overlap_full,
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro, &slicenoo,
			&rois_pixelo, &roiw_pixelo, &scalerad_pixelo, &weight_pixelo, &softrad_pixelo, &zcdelt3o,
			&spaxel_fluxo, &spaxel_weighto, &spaxel_varo, &spaxel_ifluxo, &spaxel_dqo,
			&roiw_ave,
			&cdelt1, &cdelt2)) {
    char  new_error[max_size_error] = "cube_match_sky_pointcloud: Invalid Parameters";
    strcpy(error, new_error);
    flag_error = 1; 
    goto cleanup;
  }


  // check that input parameters are valid:

  if ((cdelt1 < 0) || (cdelt2 < 0)) {
    char new_error[max_size_error] = "cdelt1 and cdelt2 must be a strictly positive number.";
    strcpy(error, new_error);
    flag_error = 1;
    goto cleanup; 
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
      (!(weight_pixel = ensure_array(weight_pixelo, &free_weight_pixel))) ||
      (!(spaxel_flux = ensure_array(spaxel_fluxo, &free_spaxel_flux))) ||
      (!(spaxel_weight = ensure_array(spaxel_weighto, &free_spaxel_weight))) ||
      (!(spaxel_var = ensure_array(spaxel_varo, &free_spaxel_var))) ||
      (!(spaxel_iflux = ensure_array(spaxel_ifluxo, &free_spaxel_flux))) ||
      (!(spaxel_dq = ensure_array_int(spaxel_dqo, &free_spaxel_dq)))){
    
    char  new_error[max_size_error] = "Failure in setting up input arrays."; 
    strcpy(error, new_error);
    flag_error = 1; 
    goto cleanup;
  }      

  npt = (int) PyArray_Size((PyObject *) coord1);
  ny = (int) PyArray_Size((PyObject *) coord2);
  nz = (int) PyArray_Size((PyObject *) wave);
  if (ny != npt || nz != npt ) {
    char  new_error[max_size_error] = "Input coordinate arrays of unequal size.";
    strcpy(error, new_error);
    flag_error = 1; 
    goto cleanup;

  }

  nxx = (int) PyArray_Size((PyObject *) xc);
  nyy = (int) PyArray_Size((PyObject *) yc);
  nwave = (int) PyArray_Size((PyObject *) zc);

  ncube = nxx * nyy * nwave;
  if (ncube ==0) {

    // 0-length input arrays. Nothing to clip. Return 0-length arrays
    char  new_error[max_size_error] = "Input Arrays have zero length.";
    strcpy(error, new_error);
    flag_error = 1; 
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
			(int *) PyArray_DATA(spaxel_dq));

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
			   (int *) PyArray_DATA(spaxel_dq));
    }
  
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
			      (double *) PyArray_DATA(spaxel_flux),
			      (double *) PyArray_DATA(spaxel_weight),
			      (double *) PyArray_DATA(spaxel_var),
			      (double *) PyArray_DATA(spaxel_iflux));

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
			     (double *) PyArray_DATA(spaxel_flux),
			     (double *) PyArray_DATA(spaxel_weight),
			     (double *) PyArray_DATA(spaxel_var),
			     (double *) PyArray_DATA(spaxel_iflux));
  }


  if (status || status1) {
    char  new_error[max_size_error] = "Error encountered in drizzling";
    strcpy(error, new_error);
    flag_error = 1; 
    goto cleanup;
  } else {
    goto cleanup;
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
  if (free_spaxel_flux) Py_XDECREF(spaxel_flux);
  if (free_spaxel_var) Py_XDECREF(spaxel_var);
  if (free_spaxel_weight) Py_XDECREF(spaxel_weight);
  if (free_spaxel_iflux) Py_XDECREF(spaxel_iflux);
  if (free_spaxel_dq) Py_XDECREF(spaxel_dq);

  if (flag_error){
    PyErr_SetString(PyExc_Exception, error);
    return NULL;    
  } else{
    return Py_BuildValue("s", "Callable C-based Cube Pointcloud");
  }
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
