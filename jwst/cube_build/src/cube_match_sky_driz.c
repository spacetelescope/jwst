/*
An IFU spectral cube is built using the overlap between the detector pixels mapped
to the 3-D spectral grid. This method is similar to 2D drizzle type mapping techniques. 
This is the top level c routine that interfaces with the python wrapper.


Main function for Python: cube_wrapper_driz

Python signature:  result = cube_wrapper_driz(
                            instrument,
                            flag_dq_plane,
                            start_region, end_region,
                            overlap_partial,
                            overlap_full,
                            xcoord, ycoord, zcoord,
                            coord1, coord2,
                            wave,
                            flux,
                            err,
                            slice_no,
                            xi1, eta1,
                            xi2, eta2,
                            xi3, eta3,
                            xi4, eta4,
                            dwave,
                            cdelt3_normal,
                            cdelt1,
                            cdelt2,
                            cdelt3_mean,
                            linear,
                            x_det,
                            y_det,
                            debug_cube_index
                        )

The output of this function is a tuple of 5 arrays:(spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq)

Parameters
----------
instrument : int
    0 = MIRI, 1 = NIRSPEC. Used for set the dq plane
flag_dq_plane : int
   0 do set the DQ plane based on FOV, but set all values =0
   1 set the DQ plane based on the FOV
start_region : int
    starting slice number for detector region used in dq flagging
end_region: int
    ending slice number for detector region used in dq flagging
overlap_partial : int
    a dq flag indicating that only a portion of the spaxel is overlapped by a mapped detector pixel
overlap_full : int
    a dq flag indicating that the entire spaxel is overlapped by the mapped detector pixel
xcoord : ndarray of float
   size of naxis1. This array holds the center x axis values of the ifu cube
ycoord : ndarray of float
   size of naxis2. This array holds the center y axis values of the ifu cube
zcoord : ndarray of float
   size of naxis3. This array holds the center z axis /wavelength values of the ifu cube
coord1 : ndarray of float
   The tangent projected xi values of the center of detector pixel.
coord2 : ndarray of float
   The tangent projected wave values of the center of detector pixel..
wave : ndarray of float
   The wavelength corresponding to the center of the detector pixel.
flux : ndarray of float
   size: point cloud elements. Flux of each point cloud member
err : ndarray of float
   size: point cloud elements. err of each point cloud member
slice_no : int
   slice number of point cloud member to be in dq flagging
xi1, eta1 : ndarray of floats
   xi, eta coordinates (tangent projection) of the detector pixel corner 1. 
xi2, eta2 : ndarray of floats
   xi, eta coordinates (tangent projection) of the detector pixel corner 2.  
xi3, eta3 : ndarray of floats
   xi, eta coordinates (tangent projection) of the detector pixel corner 3.  
xi4, eta4 : ndarray of floats
   xi, eta coordinates (tangent projection) of the detector pixel corner 4.  
dwave : ndarray of float
   Delta wavelength ranges for pixel
cdelt3_normal : ndarray of float
   For linear wavelength range = 0
   For non-linear wavelength range = delta wavelength of neighboring wavelength planes
cdelt1 : double
   Naxis 1 scale for cube
cdelt2 : double
   Naxis 2 scale for cube
cdelt3_mean: ndarray of float
   Average of cdelt3_normal
linear : bool
   Type of wavelength grid for IFU cube: linear or non-linear
x_det : ndarray of float
   size: point cloud elements. X detector value of each point cloud member
y_det : ndarray of float
   size: point cloud elements. Y detector value of each point cloud member
debug_cube_index : int
   if > 0, value of cube index to print information on 

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

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_match_sky_driz_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// routines used from cube_utils.c

extern int alloc_flux_arrays(int nelem, double **fluxv, double **weightv, double **varv,  double **ifluxv);


extern int dq_miri(int start_region, int end_region, int overlap_partial, int overlap_full,
		   int nx, int ny, int nz,
		   double cdelt1, double cdelt2, double cdelt3_mean,
		   double *xc, double *yc, double *zc,
		   double *coord1, double *coord2, double *wave,
		   double *sliceno,
		   long ncube, long npt,
		   int **spaxel_dq);

extern int dq_nirspec(int overlap_partial,
		      int nx, int ny, int nz,
		      double cdelt1, double cdelt2, double cdelt3_mean,
		      double *xc, double *yc, double *zc,
		      double *coord1, double *coord2, double *wave,
		      double *sliceno,
		      long ncube, long npt,
		      int **spaxel_dq);

extern int set_dqplane_to_zero(int ncube, int **spaxel_dq);

extern double sh_find_overlap(double xcenter, double ycenter,
                              double xlength, double ylength,
                              double xPixelCorner[],double yPixelCorner[]);

// extern double find_area_quad(double MinX, double MinY, double Xcorner[], double Ycorner[]);

// C routine that  does that does the drizzling
int match_driz(double *xc, double *yc, double *zc,
	       double *wave,
	       double *flux, double *err,
	       double *xi1, double *eta1,double *xi2, double *eta2,
	       double *xi3, double *eta3,double *xi4, double *eta4,
	       double *dwave,
	       double *cdelt3,
	       double *x_det, double *y_det,
	       double cdelt1, double cdelt2,
	       int nx, int ny, int nwave, long ncube, long npt, int linear, long debug_cube_index,
	       double **spaxel_flux, double **spaxel_weight, double **spaxel_var,
	       double **spaxel_iflux) {

  // xc : IFU grid point values along x-axis
  // yc : IFU grid point values along y-axis
  // zc : IFU grid point values along the z-axis
  // wave : wavelength of pixels
  // flux : flux values of pixels
  // err : err values of the pixels
  // xi1, eta1 :  xi, eta coordinates of a corner 1 of a pixel
  // xi2, eta2 :  xi, eta coordinates of a corner 2 of a pixel
  // xi3, eta3 :  xi, eta coordinates of a corner 3 of a pixel
  // xi4, eta4 :  xi, eta coordinates of a corner 4 of a pixel
  // dwave : delta wavelength of pixel
  // cdelt3:  wavelength sampling of IFU cube
  // x_det, y_det :x, y detector values
  // cdelt1 : spatial size along x-axis
  // cdelt2 : spatial size along y-axis
  // nx, ny, nwave : size of IFU cube in x,y,z dimensions
  // ncube : total numober cube elements
  // npt : number of detector pixels
  // linear: If T then, the IFU cube has a linear wavelength grid.
  //         If F then, the IFU cube has a non-linear wavelength grid.
  // debug_cube_index: cube index of voxel to report information on
  // spaxel_flux : return value of weighted combined flux
  // spaxel_weight : return value of combined weights
  // spaxel_var : return value of weighted combined variance
  // spaxel_iflux :return value of weighted combined iweighting map

  double *fluxv=NULL, *weightv=NULL, *varv=NULL, *ifluxv=NULL;  // vector for spaxel

  int k,j,ix1,ix2,iy1,iy2, iw1, iw2;
  int nxy, ix, iy, iw, index_xy, index_cube;
  double wdiff, zreg;
  double w1;
  double weighted_flux, weighted_var;
  double max_dwave;
  double xpixel[5], ypixel[5];
  double xmax, ymax, xmin, ymin, max_cdelt3, area, area_weight;

  double ptmin, ptmax, spxmin, spxmax, zoverlap, z1, z2, z3;
  double cdelt1_half, cdelt2_half;
  double xleft, xright, ybot, ytop;
  // double area_quad;
  // allocate memory to hold output
  if (alloc_flux_arrays(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;


  // find max of cdelt3, dwave to be used to estimate which wavelength plane the
  // pixel falls on
  zreg =0;
  max_cdelt3 = cdelt3[0];
  max_dwave = dwave[0];
  for (iw =1; iw < nwave;  iw++){
    if(cdelt3[iw] > max_cdelt3) { max_cdelt3 = cdelt3[iw];}
    if(dwave[iw] > max_dwave){ max_dwave = dwave[iw];}
  }

  // printf("debug_spaxel  %i  \n ", debug_cube_index);
  // loop over each detector pixel and find which spaxels it overlaps with
  nxy = nx * ny;
  for (k = 0; k < npt; k++) {
      xpixel[0] = xi1[k];
      xpixel[1] = xi2[k];
      xpixel[2] = xi3[k];
      xpixel[3] = xi4[k];
      xpixel[4] = xi1[k];

      ypixel[0] = eta1[k];
      ypixel[1] = eta2[k];
      ypixel[2] = eta3[k];
      ypixel[3] = eta4[k];
      ypixel[4] = eta1[k];

      xmin = xpixel[0];
      ymin = ypixel[0];
      xmax = xmin;
      ymax = ymin;
      for (j =1; j< 5; j++){
	if(xpixel[j] > xmax) xmax = xpixel[j];
	if(ypixel[j] > ymax) ymax = ypixel[j];
	if(xpixel[j] < xmin) xmin = xpixel[j];
	if(ypixel[j] < ymin) ymin = ypixel[j];
      }

      cdelt1_half = cdelt1/2.0;
      cdelt2_half = cdelt2/2.0;

      // find the area of the pixel (quadrilateral) not needed now - keeping if needed later
      // area_quad = find_area_quad(xmin, ymin, xpixel, ypixel);

      // convert to integer values to get the approximate region to search
      // cdelt1_half and cdelt2_half - may not be needed.
      ix1 = fabs((xmin - cdelt1_half - xc[0])/cdelt1)-1;
      ix2 = fabs((xmax + cdelt1_half - xc[0])/cdelt1)+1;

      iy1 = fabs((ymin - cdelt2_half - yc[0])/cdelt2)-1;
      iy2 = fabs((ymax + cdelt2_half - yc[0])/cdelt2)+1;
      if(ix1 < 0) { ix1 = 0;}
      if(iy1 < 0) { iy1 = 0;}
      if(ix2 > nx) { ix2 = nx;}
      if(iy2 > ny) { iy2 = ny;}

      // estimate the wavelength overlapping region using max_cdelt3 and max_dwave
      // estimating wavelength range works if we have a linear wavelength
      if (linear == 1){
	w1 = wave[k] - (max_cdelt3 + max_dwave) - zc[0];
	if(w1 < 0){
	  iw1 = 0;
	}else{
	  iw1 = fabs((w1)/ (max_cdelt3+max_dwave));
	}
	iw2 = ceil(fabs((wave[k] + (max_cdelt3 + max_dwave) - zc[0])/max_cdelt3));
      } else{
	iw1 = 0;
	iw2 = nwave;
      }

      iw1 = 0;
      iw2 = nwave;
      for (iw =iw1; iw < iw2;  iw++){
	zreg = fabs(dwave[k] + cdelt3[iw]);
	// zreg = 0.0025; (roiw size for testing- usually larger than zreg)
	wdiff = zc[iw] - wave[k];

	if( fabs(wdiff) < zreg){
	  // Fractional wavelength overlaps to use for weighting
	  ptmin = wave[k] - dwave[k]/2;
	  ptmax = wave[k] + dwave[k]/2;
	  spxmin = zc[iw] - cdelt3[iw]/2;
	  spxmax = zc[iw] + cdelt3[iw]/2;
	  z1 = spxmax - ptmin;
	  z2 = spxmax - ptmax;
	  z3 = spxmin - ptmin;
	  if(z1 < 0) { z1 =0;}
	  if(z2 < 0) { z2 =0;}
	  if(z3 < 0) { z3 =0;}
	  zoverlap = z1 - z2 - z3;
	  if(zoverlap < 0) { zoverlap = 0;}

	  // find match in spatial dimension using approximate locations based on
	  // ix1, ix2, iy1, iy2
	  for (ix =ix1; ix < ix2; ix++){
	    for (iy =iy1; iy < iy2; iy++){

	      // narrow down the spatial region
	      xleft = xc[ix] - cdelt1*0.5;
	      xright = xc[ix] + cdelt1*0.5;

	      ybot = yc[iy] - cdelt2*0.5;
	      ytop = yc[iy] + cdelt2*0.5;

	      index_xy = iy* nx + ix;

	      if(xleft < xmax && xright > xmin && ybot < ymax && ytop > ymin){

		index_xy = iy* nx + ix;
		index_cube = iw*nxy + index_xy;
		// Spatial overlap between detector pixel and cube spaxel
		area = sh_find_overlap(xc[ix], yc[iy],
				       cdelt1, cdelt2,
				       xpixel,ypixel);

		// area_weight = area of overlap * wavelength overlap
		area_weight = area * zoverlap;
		if(area_weight > 0) {
		  weighted_flux =  flux[k]* area_weight;
		  weighted_var = (err[k]* area_weight) * (err[k]*area_weight);
		  fluxv[index_cube] = fluxv[index_cube] + weighted_flux;
		  weightv[index_cube] = weightv[index_cube] + area_weight;
		  varv[index_cube] = varv[index_cube] + weighted_var;
		  ifluxv[index_cube] = ifluxv[index_cube] +1.0;
		}
		
		// Keep print statement in code - used for debugging
		if (index_cube == debug_cube_index){
		  printf("spaxel, flux, x, y [count starting at 0]  %i %f %f %f  \n ", index_cube,
			 flux[k], x_det[k], y_det[k]);
		}
		
		
		// end of print statements
	      } // xleft, xright, ybot, ytop

	    }// end loop over iy
	  } // end loop over ix
	} // check of wave
      } // end loop over wave
  } // end loop over detector elements

    // assign output values:
  *spaxel_flux = fluxv;
  *spaxel_weight = weightv;
  *spaxel_var = varv;
  *spaxel_iflux = ifluxv;

  return 0;
}


// This is a C structure that represents a NumPy array in the C API.
// This ensures that all the numpy arrays passed to the C routines
// follow C array rules. 
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


// Wrapper code that is called from python code and sets up interface with C code.
// cube_wrapper_driz.

static PyObject *cube_wrapper_driz(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xco, *yco, *zco, *fluxo, *erro, *coord1o, *coord2o, *waveo, *slicenoo; // codespell:ignore erro
  PyObject  *cdelt3o;
  PyObject *xi1o, *eta1o, *xi2o, *eta2o, *xi3o, *eta3o, *xi4o, *eta4o, *x_deto, *y_deto;
  PyObject *dwaveo;
  
  double cdelt1, cdelt2,cdelt3_mean;
  int  nwave, nxx, nyy;
  long npt, ncube, debug_cube_index;
  int linear;
  int instrument, flag_dq_plane, start_region, end_region, overlap_partial, overlap_full;
  double *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  double *spaxel_iflux=NULL;
  int *spaxel_dq=NULL;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_flux=0, free_err=0, free_cdelt3=0;
  int free_sliceno=0, free_x_det=0, free_y_det=0;
  int free_xi1=0, free_eta1=0,free_xi2=0, free_eta2=0,free_xi3=0, free_eta3=0,free_xi4=0, free_eta4=0;
  int free_dwave=0;

  PyArrayObject *xc, *yc, *zc, *flux, *err, *coord1, *coord2, *wave;
  PyArrayObject *xi1, *eta1, *xi2, *eta2, *xi3, *eta3, *xi4, *eta4, *dwave;
  PyArrayObject *cdelt3, *sliceno, *x_det, *y_det;
  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_var_arr=NULL;
  PyArrayObject *spaxel_iflux_arr=NULL, *spaxel_dq_arr=NULL;
  npy_intp npy_ncube = 0;

  int  ny,nz;

  if (!PyArg_ParseTuple(args, "iiiiiiOOOOOOOOOOOOOOOOOOOdddiOOl:cube_wrapper_driz",
			&instrument, &flag_dq_plane, &start_region, &end_region, &overlap_partial, &overlap_full,
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro, &slicenoo, // codespell:ignore erro
			&xi1o, &eta1o, &xi2o, &eta2o, &xi3o, &eta3o, &xi4o, &eta4o,
			&dwaveo,
			&cdelt3o, &cdelt1, &cdelt2, &cdelt3_mean, &linear,
			&x_deto, &y_deto, &debug_cube_index)) {
    return NULL;
  }


  // check that input parameters are valid:

  if ((cdelt1 <= 0) || (cdelt2 <= 0)) {
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
      (!(x_det = ensure_array(x_deto, &free_x_det))) ||
      (!(y_det = ensure_array(y_deto, &free_y_det))) ||
      (!(err = ensure_array(erro, &free_err))) || // codespell:ignore erro
      (!(sliceno = ensure_array(slicenoo, &free_sliceno))) ||
      (!(xi1 = ensure_array(xi1o, &free_xi1))) ||
      (!(eta1 = ensure_array(eta1o, &free_eta1))) ||
      (!(xi2 = ensure_array(xi2o, &free_xi2))) ||
      (!(eta2 = ensure_array(eta2o, &free_eta2))) ||
      (!(xi3 = ensure_array(xi3o, &free_xi3))) ||
      (!(eta3 = ensure_array(eta3o, &free_eta3))) ||
      (!(xi4 = ensure_array(xi4o, &free_xi4))) ||
      (!(eta4 = ensure_array(eta4o, &free_eta4))) ||
      (!(cdelt3 = ensure_array(cdelt3o, &free_cdelt3))) ||
      (!(dwave = ensure_array(dwaveo, &free_dwave))) )
    {
      goto cleanup;

    }

  npt = (long) PyArray_Size((PyObject *) coord1);
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

  // if flag_dq_plane = 1, Set up the dq plane of the IFU cube. 
  
  int status1 = 0;

  if(flag_dq_plane){
    if (instrument == 0){
      status1 = dq_miri(start_region, end_region,overlap_partial, overlap_full,
			nxx,nyy,nwave,
			cdelt1, cdelt2, cdelt3_mean,
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
			   cdelt1, cdelt2, cdelt3_mean,
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
  // Driz the mapped detector data onto the IFU cube
  //______________________________________________________________________
    status = match_driz((double *) PyArray_DATA(xc),
			(double *) PyArray_DATA(yc),
			(double *) PyArray_DATA(zc),
			(double *) PyArray_DATA(wave),
			(double *) PyArray_DATA(flux),
			(double *) PyArray_DATA(err),
			(double *) PyArray_DATA(xi1),
			(double *) PyArray_DATA(eta1),
			(double *) PyArray_DATA(xi2),
			(double *) PyArray_DATA(eta2),
			(double *) PyArray_DATA(xi3),
			(double *) PyArray_DATA(eta3),
			(double *) PyArray_DATA(xi4),
			(double *) PyArray_DATA(eta4),
			(double *) PyArray_DATA(dwave),
			(double *) PyArray_DATA(cdelt3),
			(double *) PyArray_DATA(x_det),
			(double *) PyArray_DATA(y_det),
			cdelt1, cdelt2,
			nxx, nyy, nwave, ncube, npt,linear, debug_cube_index,
			&spaxel_flux, &spaxel_weight, &spaxel_var, &spaxel_iflux);


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
  if (free_cdelt3) Py_XDECREF(cdelt3);
  if (free_sliceno) Py_XDECREF(sliceno);
  if (free_xi1) Py_XDECREF(xi1);
  if (free_xi2) Py_XDECREF(xi2);
  if (free_xi3) Py_XDECREF(xi3);
  if (free_xi4) Py_XDECREF(xi4);
  if (free_eta1) Py_XDECREF(eta1);
  if (free_eta2) Py_XDECREF(eta2);
  if (free_eta3) Py_XDECREF(eta3);
  if (free_eta4) Py_XDECREF(eta4);
  if (free_dwave) Py_XDECREF(dwave);
  return result;
}


static PyMethodDef cube_methods[] =
{
    {
        "cube_wrapper_driz",
	cube_wrapper_driz,
	METH_VARARGS,
        "cube_wrapper(put in doc string)"
    },
    {NULL, NULL, 0, NULL}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cube_match_sky_driz",             /* m_name */
    "find point cloud matches for each spaxel center",  /* m_doc */
    -1,                          /* m_size */
    cube_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_cube_match_sky_driz(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
