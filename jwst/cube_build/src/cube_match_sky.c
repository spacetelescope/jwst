/*
The detector pixels are represented by a 'point could' on the sky. The IFU cube is
represented by a 3-D regular grid. This module finds the point cloud members contained
in a region centered on the center of the cube spaxel. The size of the spaxel is spatial
coordinates is cdetl1 and cdelt2, while the wavelength size is zcdelt3.
This module uses the e modified shephard weighting method to determine how to  weight each point clold member
in the spaxel.
 
Main function for Python: cube_wrapper

Python signature: result = cube_wrapper(instrument, start_region, end_region,
                                        self.overlap_partial, self.overlap_full,
                                        self.xcoord, self.ycoord, self.zcoord,
                                        coord1, coord2, wave, flux, err, slice_no,
                                        rois_pixel, roiw_pixel, scalerad_pixel,self.cdelt3_normal,
                                        roiw_ave, self.cdelt1, self.cdelt2)
provide more details

The output of this function is a tuple of 5 arrays:(spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq) 
example output 

Parameters
----------
instrument: int
    0 = MIRI, 1 = NIRSPEC. Used for set the dq plane
start_region : int
    starting slice number for detector region
end_region: int 
    ending slice number for detector region
xcoord : numpy.ndarray
   size of naxis1. This array holds the center x axis values of the ifu cube 
ycoord : numpy.ndarray
   size of naxis2. This array holds the center y axis values of the ifu cube 
zcoord : numpy.ndarray
   size of naxis3. This array holds the center x axis values of the ifu cube 
flux : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
err : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
coord1 : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
coord2 : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
wave : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
rois_pixel : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
roiw_pixel : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
scalerad_pixel : numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber
zcdelt3: numpy.ndarray
   size: point cloud elements. Flux of each point cloud memeber

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
spaxel_dq : numpy.ndarray

*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
//#include "cube_utils.c"

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_cube_match_sky_numpy_api    //WHAT IS THIS AND WHERE IS IT USED???
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION



// routines used from cube_utils.c

extern double sh_find_overlap(const double xcenter, const double ycenter, 
		      const double xlength, const double ylength,
			      double xPixelCorner[],double yPixelCorner[]);

// Allocate the memory to for the spaxel arrays in the cube

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

//________________________________________________________________________________
// allocate the memory for the spaxel DQ array

int mem_alloc_dq(long nelem, int **idqv) {
    int  *i;
    const char *msg = "Couldn't allocate memory for output arrays.";

    i = (int*)malloc(nelem * sizeof(int));
    if (i) {
        *idqv = i;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }
    return 0;
}


//________________________________________________________________________________
// Routine for MIRI DQ plane assignment

int corner_wave_plane_miri(int w, int start_region, int end_region,
			   double roiw_ave,
			   double *zc,
			   double *coord1, double *coord2, double *wave,
			   double *sliceno,
			   long ncube, int npt,
			   double *corner1, double *corner2, double *corner3, double *corner4) {
  /* 
     For wavelength plane determine the corners (in xi,eta) of the FOV for MIRI
     w : wavelength plane
     start_region : starting slice # for channel (slice # in nirspec) 
     end_region : ending slice # for chanel      (slice # in nirspec)
     roiw_ave : average roiw for all wavelengths
     zc : array of wavelenths
     coord1 : point cloud xi values
     coord2 : point cloud eta values
     wave : point cloud wavelength values
     sliceno : point cloud slice no
     ncube : number of cube values
     npt: number of point cloud elements
     corner1 : xi, eta of corner 1
     corner2 : xi, eta of corner 2
     corner3 : xi, eta of corner 3
     corner4 : xi, eta of corner 4
   */

  int status = 0; 
  int ic1_start_min = -1;
  int ic1_start_max = -1;
  int ic2_start_min = -1;
  int ic2_start_max = -1;

  float c1_start_min = 10000;
  float c2_start_min = 10000;
  float c1_start_max = -10000;
  float c2_start_max = -10000;

  int ic1_end_min = -1;
  int ic1_end_max = -1;
  int ic2_end_min = -1;
  int ic2_end_max = -1;

  float c1_end_min = 10000;
  float c2_end_min = 10000;
  float c1_end_max = -10000;
  float c2_end_max = -10000;

  for (int ipt =0; ipt< npt ; ipt++){
    int slice = (int)sliceno[ipt];
    double wave_distance = fabs(zc[w] - wave[ipt]);


    float c11 = -1;  // coord1 for start region
    float c21 = -1;  // coord2 for start region
    float c12 = -1;  // coord1 for end region
    float c22 = -1;  // corrd2 for end region

    if(slice == start_region || slice==end_region){
      // Find all the coordinates on wave slice with slice = start region
      // These points will define corner 1 (min c2) and corner 2 (max c2)
      if (wave_distance < roiw_ave && slice == start_region){ 
	c11 = coord1[ipt];
	c21 = coord2[ipt];
      }

      // Find all the coordinates on wave slice with slice = start region
      // These points will define corner 4 (min c2) and corner 3 (max c2)
      if (wave_distance < roiw_ave && slice == end_region){
	c12 = coord1[ipt];
	c22 = coord2[ipt];
      }

      // for the start region find min and max c1,c2
      if( c11 !=-1){
	// check c1 points for min
	if (c11 < c1_start_min) {
	  c1_start_min = c11;
	  ic1_start_min = ipt;
	}
	// check c1 points for max
	if (c11 > c1_start_max) {
	  c1_start_max= c11;
	  ic1_start_max = ipt;
	}
	// check c2 points for min
	if (c21 < c2_start_min) {
	  c2_start_min = c21;
	  ic2_start_min = ipt;
	}
	// check c2 points for max
	if (c21 > c2_start_max) {
	  c2_start_max= c21;
	  ic2_start_max = ipt;
	}
      }
      // for the end region find min and max c1,c2
      if( c12 !=-1){
	if (c12 < c1_end_min) {
	  c1_end_min = c12;
	  ic1_end_min = ipt;
	}
	if (c12 > c1_end_max) {
	  c1_end_max = c12;
	  ic1_end_max = ipt;
	}
	if (c22 < c2_end_min) {
	  c2_end_min = c22;
	  ic2_end_min = ipt;
	}
	if (c22 > c2_end_max) {
	  c2_end_max = c22;
	  ic2_end_max = ipt;
	}
      }
    }
  } // end looping over point cloud

  // Make sure the 2 extreme slices are found on the FOV. Not finding both can occur for edge wavelength planes
  // or empty wavelength planes between channels

  if( ic1_start_min ==-1 || ic1_start_max ==-1 || ic1_end_min == -1 || ic1_end_max == -1){
    status = 1;
    return status;
  } else {
  // Find the length in the start and end regions for c1 and c2. This will help define which coordinates to use to set corners.
  // Because we do not know the orientation on the sky pick the  longest length to set how to pick corners. 
    float length_c1_start = c1_start_max - c1_start_min;
    float length_c2_start = c2_start_max - c2_start_min;
    float length_c1_end = c1_end_max - c1_end_min;
    float length_c2_end = c2_end_max - c2_end_min;

    int c1_use = 1; // use the c1 coords to set corners 
    if(length_c1_start < length_c2_start){
      c1_use = 0;   // use the c2 coords to set the corners
    }
    if(c1_use ==0) {
      if (length_c1_end > length_c2_end){
	printf(" Check sizes 1 \n" );
	printf(" for wavelength %i start  %f %f %f %f  \n", w,c1_start_min, c1_start_max, c2_start_min, c2_start_max);
	printf(" for wavelength %i start %f %f %f %f  \n", w,c1_end_min, c1_end_max, c2_end_min, c2_end_max);
	printf("length %f %f %f %f \n ", length_c1_start, length_c2_start, length_c1_end, length_c2_end);
      }
    }
      
    if(c1_use == 1) {
      if (length_c1_end < length_c2_end){
	printf(" Check sizes 2 \n " );
	printf(" for wavelength %i start  %f %f %f %f  \n", w,c1_start_min, c1_start_max, c2_start_min, c2_start_max);
	printf(" for wavelength %i start %f %f %f %f  \n", w,c1_end_min, c1_end_max, c2_end_min, c2_end_max);
	printf("length %f %f %f %f \n ", length_c1_start, length_c2_start, length_c1_end, length_c2_end);
      }
    }

    if(c1_use ==0) {
      corner1[0] = coord1[ic2_start_min];
      corner1[1] = coord2[ic2_start_min];
    
      corner2[0] = coord1[ic2_start_max];
      corner2[1] = coord2[ic2_start_max];
    
      corner3[0] = coord1[ic2_end_max];
      corner3[1] = coord2[ic2_end_max];

      corner4[0] = coord1[ic2_end_min];
      corner4[1] = coord2[ic2_end_min];
    } else{
      corner1[0] = coord1[ic1_start_min];
      corner1[1] = coord2[ic1_start_min];
    
      corner2[0] = coord1[ic1_start_max];
      corner2[1] = coord2[ic1_start_max];
    
      corner3[0] = coord1[ic1_end_max];
      corner3[1] = coord2[ic1_end_max];
    
      corner4[0] = coord1[ic1_end_min];
      corner4[1] = coord2[ic1_end_min];
    }

    status = 0;
    return status;
  }
}

//________________________________________________________________________________
// MIRI DQ routine. Find the overlap of the FOV for the wavelength slice with sky
int overlap_fov_with_spaxels(int overlap_partial,  int overlap_full,
                             double cdelt1, double cdelt2,
                             int naxis1, int naxis2,
                             double xcenters[], double ycenters[],
                             double xi_corner[], double eta_corner[],
			     int wave_slice_dq[]) {

  /* find the amount of overlap of FOV  each spaxel

        Given the corners of the FOV  find the spaxels that
        overlap with this FOV.  Set the intermediate spaxel  to
        a value based on the overlap between the FOV for each exposure
        and the spaxel area. The values assigned are:
        a. overlap_partial = overlap partial
        b  overlap_full = overlap_full
        bit_wise combination of these values is allowed to account for
        dithered FOVs.

        Parameter
        ----------
        xi_corner: xi coordinates of the 4 corners of the FOV on the wavelenghth plane
        eta_corner: eta coordinates of the 4 corners of the FOV on the wavelength plane
        wmin: minimum wavelength bin in the IFU cube that this data covers
        wmax: maximum wavelength bin in the IFU cube that this data covers

        Sets
        -------
        wave_slice_dq: array containing intermediate dq flag

  */

  // loop over spaxels in the wavelength plane and set slice_dq
  // roughly find the spaxels that might be overlapped

  int status = 0;
  double ximin = 1000.0;
  double etamin = 1000.0;
  double ximax = -10000.0;
  double etamax = -1000.0;
  for (int i=0; i<4 ; i++){
    if (xi_corner[i] < ximin){ ximin = xi_corner[i];}
    if (xi_corner[i] > ximax){ ximax = xi_corner[i];}
    if (eta_corner[i] < etamin){ etamin = eta_corner[i];}
    if (eta_corner[i] < etamax){ etamax = eta_corner[i];}
  }

  double area_box = cdelt1 * cdelt2;
  float tolerance_dq_overlap = 0.05;  //spaxel has to have 5% overlap to flag in FOV
  // loop over cube xcenters and cube ycenters
  for (int ix = 0; ix < naxis1; ix++){
    double x1 = (xcenters[ix] - cdelt1)/2;
    double x2 = (xcenters[ix] + cdelt1)/2;
    if(x1 > ximin && x2 < ximax){
      for (int iy = 0; iy< naxis2; iy++){

	double y1 = (ycenters[ix] - cdelt2)/2;
	double y2 = (ycenters[ix] + cdelt2)/2;
	if(y1 > etamin && y2 < etamax){
	  int ixy = iy* naxis1 + ix;
	  double area_overlap = sh_find_overlap(xcenters[ix],ycenters[iy],
						cdelt1, cdelt2,
						xi_corner, eta_corner);

	  double overlap_coverage = area_overlap / area_box;
	  
	  if (overlap_coverage > tolerance_dq_overlap){
	    if (overlap_coverage > 0.95) {
	      wave_slice_dq[ixy] = overlap_full;
	    }else{
	      wave_slice_dq[ixy] = overlap_partial;
	    }
	  }
	}// end y1, y2 test
      }// end loop over iy
    } // end x1, x2 test
  } // end loop over ix
  
  return status ;
 
}


//________________________________________________________________________________
// Routine to setting NIRSpec dq plane for each wavelength plane

int slice_wave_plane_nirspec(int w, int slicevalue,
		      double roiw_ave,
		      double *zc,
		      double *coord1, double *coord2, double *wave,
		      double *sliceno,
		      long ncube, int npt,
		      double *c1_min, double *c1_max, double *c2_min, double *c2_max) {
  /* 
     For wavelength plane determine the limits of each slice
     w : wavelength plane
     slicevalue : slice # 1 to 30 
     roiw_ave : average roiw for all wavelengths
     zc : array of wavelenths
     coord1 : point cloud xi values
     coord2 : point cloud eta values
     wave : point cloud wavelength values
     sliceno : point cloud slice no - starts at 1
     ncube : number of cube values
     npt: number of point cloud elements

     return:
     c1_min, c1_max, c2_min, c2_max
   */

  double dvalue = 10000;
  *c1_min = dvalue;
  *c2_min = dvalue;
  *c1_max = -dvalue;
  *c2_max = -dvalue;
 
  for (int ipt =0; ipt< npt ; ipt++){
    int slice = (int)sliceno[ipt];
    double wave_distance = fabs(zc[w] - wave[ipt]);

    // Find all the coordinates on wave slice with slice = start region

    if (wave_distance < roiw_ave && slice == slicevalue){

      double c1 = coord1[ipt];
      double c2 = coord2[ipt];
      // find min, max of xi eta
      if (c1 < *c1_min) {*c1_min = c1;}
      if (c1 > *c1_max) {*c1_max = c1;}
      if (c2 < *c2_min) {*c2_min = c2;}
      if (c2 > *c2_max) {*c2_max = c2;}
    }

  } // end looping over point cloud

  int status = 0;
  if(*c1_min == dvalue || *c2_min == dvalue || *c1_max ==-dvalue || *c2_max == -dvalue){
    // Problem finding limits of slice for wavelength plane
    // This is likely caused the no valid data on wavelength plane
    // The two ends of wavelengths can have DQ detector data set to DO_NOT_USE - setting up no
    // valid data on thee planes in the IFU Cube. 

    status = 1; 
  }

  //if (status){
  //  printf(" found min max of slice for wavelength %f %i %f %f %f %f \n" , zc[w],slicevalue, *c1_min, *c1_max, *c2_min, *c2_max );
  //}
  return status;
}



//________________________________________________________________________________

// NIRSpec  Find the overlap of the slices  for the wavelength slice with sky
int overlap_slice_with_spaxels(int overlap_partial,
			       double cdelt1, double cdelt2,
			       int naxis1, int naxis2,
			       double xcenters[], double ycenters[],
			       double xi_min, double eta_min,
			       double xi_max, double eta_max,
			       int wave_slice_dq[]) {

  /* 
     Set the initial dq plane indicating if the input data falls on a spaxel

     This algorithm assumes the input data falls on a line in the IFU cube, which is
     the case for NIRSpec slices. The NIRSpec slice's endpoints are used to determine
     which IFU spaxels the slice falls on to set an initial dq flag.
     
     Parameters
     ---------
     overlap_partial: intermediate dq flag

     wavelength: the wavelength bin of the IFU cube working with

     Sets
     ----
     wave_slice_dq : numpy.ndarray containing intermediate dq flag

     Bresenham's Line Algorithm to find points a line intersects with grid.

     Given the endpoints of a line find the spaxels this line intersects.

     Returns
     -------
     Points: a tuple of x,y spaxel values that this line intersects

  */

  int status = 0;

  //set up line - convert to integer values
  int x1 = (int)((xi_min - xcenters[0]) / cdelt1);
  int y1 = (int)((eta_min - ycenters[0]) / cdelt2);
  int x2 = (int)((xi_max - xcenters[0]) / cdelt1);
  int y2 = (int)((eta_max - ycenters[0]) / cdelt2);

  int dx = x2 - x1;
  int dy = y2 - y1;

  bool is_steep;
  is_steep = abs(dy) > abs(dx);

  // if is_steep switch x and y 
  if (is_steep){
      x1 = y1;
      y1 = x1;
      x2 = y2;
      y2 = x2;
    }

  // Swap start and end points if necessary and store swap state
  bool swapped;
  swapped = false;
  if (x1 > x2){
    x1 = x2;
    x2 = x1;

    y1 = y2;
    y2 = y1;
    swapped = true;
  }

  // Recalculate differences
  dx = x2 - x1;
  dy = y2 - y1;

  //calculate error
  int error = (int)(dx / 2.0);
  int ystep = -1;
  if (y1 < y2){
    ystep = 1;
  }

  //printf(" corners %f %f %f %f \n ", xi_min, xi_max, eta_min, eta_max);
  //printf(" for slice on wavelength x y ranges on cube %i %i %i %i \n ", x1, x2, y1, y2);
  // iterate over grid to generate points between the start and end of line
  int y = y1;
  for (int x = x1; x< (x2 + 1); x++){
    int yuse  = y;
    int xuse = x ;
    if (is_steep){
	yuse = x;
	xuse = y;
      }
    //coord = (y, x) if is_steep else (x, y);
    int index = (yuse * naxis1) + xuse;
    wave_slice_dq[index] = overlap_partial;
    error -= abs(dy);
    if (error < 0){
      y += ystep;
      error += dx;
    }

  }
  return status; 
}

//________________________________________________________________________________
int dq_miri(int start_region, int end_region, int overlap_partial, int overlap_full,
	    int nx, int ny, int nz,
	    double cdelt1, double cdelt2, double roiw_ave,
	    double *xc, double *yc, double *zc,
	    double *coord1, double *coord2, double *wave,
	    double *sliceno,
	    long ncube, int npt, 
	    int **spaxel_dq) {

  int *idqv = NULL;  // int vector for spaxel

    if (mem_alloc_dq(ncube, &idqv)) return 1;

    // Set all data to zero
    for (long i = 0; i < ncube; i++){
      idqv[i] = 0;
    }

    double corner1[2];
    double corner2[2];
    double corner3[2];
    double corner4[2];
    
    // for each wavelength plane find the 2 extreme slices to set FOV. Use these two extreme slices to set up the
    // corner of the FOV for each wavelength

    int nxy = nx * ny;
    // Loop over the wavelength planes and set DQ plane 
    for (int w = 0; w  < nz; w++) {

      int wave_slice_dq[nxy];
      for( int i = 0; i < nxy; i ++){
	wave_slice_dq[i] = 0;
      }
      
      int status =  corner_wave_plane_miri( w, start_region, end_region, roiw_ave, zc,
				   coord1, coord2, wave, sliceno, ncube, npt,
				   corner1, corner2, corner3, corner4);
      if( status == 0){ // found min and max slice on wavelengh plane
	double xi_corner[4];
	double eta_corner[4];
	xi_corner[0] = corner1[0];
	xi_corner[1] = corner2[0];
	xi_corner[2] = corner3[0];
	xi_corner[3] = corner4[0];
	
	eta_corner[0] = corner1[1];
	eta_corner[1] = corner2[1];
	eta_corner[2] = corner3[1];
	eta_corner[3] = corner4[1];
	


	status = overlap_fov_with_spaxels(overlap_partial, overlap_full,
					  cdelt1,cdelt2,
					  nx, ny,
					  xc, yc,
					  xi_corner, eta_corner,
					  wave_slice_dq);
	long istart = nxy*w;
	long iend = istart + nxy;
	for( long in = istart; in < iend; in ++){
	  long ii = in - istart;
	  idqv[in] = wave_slice_dq[ii];
	}
      }
    } // end loop over wavelength

    *spaxel_dq = idqv;

    return 0;
}


//________________________________________________________________________________
// set up the dq plane for NIRSPEC
 
int dq_nirspec(int overlap_partial,
	       int nx, int ny, int nz,
	       double cdelt1, double cdelt2, double roiw_ave,
	       double *xc, double *yc, double *zc,
	       double *coord1, double *coord2, double *wave,
	       double *sliceno,
	       long ncube, int npt,
	       int **spaxel_dq) {

  /*
    Set an initial DQ flag for the NIRSPEC IFU cube based on FOV of input data.

    Map the FOV of each NIRSpec slice  to the DQ plane and set an initial DQ
    flagging. For  NIRSpec the 30 different slices map to different FOV on the
    range of wavelengths.  The FOV of the slice is really just a line, so instead of using
    the routines the finds the overlap between a polygon and regular grid-
    which is used for MIRI - an algorithm that determines the spaxels that
    the slice line intersects is used instead.
    
    Paramteter
    ---------
    overlap_partial - intermediate DQ flag
    coord1: xi coordinates of input data (~x coordinate in IFU space)
    coord2: eta coordinates of input data (~y coordinate in IFU space)
    wave: wavelength of input data
    roiw_ave: average spectral roi used to determine which wavelength bins
    the input values would be mapped to
    slice_no: integer slice value of input data (used in MIRI case to find
    the points of the edge slices.)

  */
  

  int *idqv = NULL;  // int vector for spaxel

  if (mem_alloc_dq(ncube, &idqv)) return 1;

  // Set all data to zero
  //printf( " number of elements allocating memory for %lu \n" , ncube);
  for (long i = 0; i < ncube; i++){
    idqv[i] = 0;
  }

  //  for each of the 30 slices - find the projection of this slice
  //     onto each of the IFU wavelength planes.


  for (int w = 0; w  < nz; w++) {
    int n_slice_found = 0;
    for (int islice =1; islice< 31 ; islice++){
      double c1_min;
      double c2_min;
      double c1_max;
      double c2_max;
      int status = 0; 
      status =  slice_wave_plane_nirspec( w, islice, roiw_ave, zc,
				   coord1, coord2, wave, sliceno, ncube, npt,
				   &c1_min, &c2_min, &c1_max, &c2_max);
	
      if( status ==0){
	
	n_slice_found++;
	int nxy = nx * ny;
	int wave_slice_dq[nxy];
	for (int j =0; j< nxy; j++){
	  wave_slice_dq[j] = 0;
	}
	status = overlap_slice_with_spaxels(overlap_partial,
					    cdelt1,cdelt2,
					    nx, ny,
					    xc, yc,
					    c1_min, c2_min, c1_max, c2_max,
					    wave_slice_dq);
      
	long istart = nxy*w;
	long iend = istart + nxy;
	for( long in = istart; in < iend; in ++){
	  long ii = in - istart;
	  idqv[in] = wave_slice_dq[ii];
	}
      } // end loop over status
      
    } // end loop over slices
    //if (n_slice_found ==0){
    //  printf(" no slices found on wavelength plane %i %f \n",w, zc[w]);
    //}
	
  } // end of wavelength
  *spaxel_dq = idqv;

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
    if (mem_alloc(ncube, &fluxv, &weightv, &varv, &ifluxv)) return 1;
    
    double set_zero=0.0;
    // Set all data to zero
    for (int i = 0; i < ncube; i++){
      varv[i] = set_zero;
      fluxv[i] = set_zero;
      ifluxv[i] = set_zero;
      weightv[i] = set_zero;
    }
    
    // loop over each point cloud member and find which roi spaxels it is found

    for (int k = 0; k < npt; k++) {
       // Search wave and find match
      int iwstart = -1;
      int iwend = -1;
      int ii = 0;
      int done_search_w = 0;

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
		//double wdist = fabs( wave[k] - zc[iw]);
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

		//if( ix == 24 && iy == 24 && iw ==0) {
		//  printf("found element %i %i %f %f %f  %f %f \n", k, index_cube, err[k], flux[k], ifluxv[index_cube], varv[index_cube], fluxv[index_cube]);
		//}
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
  PyObject *rois_pixelo, *roiw_pixelo, *scalerad_pixelo, *zcdelt3o;
  
  double cdelt1, cdelt2, roiw_ave;
  int  nwave, npt, nxx, nyy, ncube;

  int instrument, start_region, end_region, overlap_partial, overlap_full;
  double *spaxel_flux=NULL, *spaxel_weight=NULL, *spaxel_var=NULL;
  double *spaxel_iflux=NULL;
  int *spaxel_dq=NULL;

  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0 , free_wave=0, status=0;
  int free_rois_pixel=0, free_roiw_pixel=0, free_scalerad_pixel=0, free_flux=0, free_err=0, free_zcdelt3=0;
  int free_sliceno=0;
  
  PyArrayObject *xc, *yc, *zc, *flux, *err, *coord1, *coord2, *wave, *rois_pixel, *roiw_pixel, *scalerad_pixel;
  PyArrayObject *zcdelt3, *sliceno;
  PyArrayObject *spaxel_flux_arr=NULL, *spaxel_weight_arr=NULL, *spaxel_var_arr=NULL;
  PyArrayObject *spaxel_iflux_arr=NULL, *spaxel_dq_arr=NULL; 
  npy_intp npy_ncube = 0;

  int  ny,nz;

  if (!PyArg_ParseTuple(args, "iiiiiOOOOOOOOOOOOOddd:cube_wrapper",
			&instrument, &start_region, &end_region, &overlap_partial, &overlap_full,
			&xco, &yco, &zco, &coord1o, &coord2o, &waveo,  &fluxo, &erro, &slicenoo,
			&rois_pixelo, &roiw_pixelo, &scalerad_pixelo,&zcdelt3o, &roiw_ave,
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
      (!(scalerad_pixel = ensure_array(scalerad_pixelo, &free_scalerad_pixel))) )
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
  
  //printf(" sizes %i %i %i %i %i \n ", nxx, nyy, nwave, npt, ncube);

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

    result = Py_BuildValue("(OOOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
			   spaxel_iflux_arr, spaxel_dq_arr);

    goto cleanup;

  }

  //______________________________________________________________________
  // First set up the dq plane
  //______________________________________________________________________
  int status1 = 0; 
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

  //______________________________________________________________________
  // Match the point cloud elements to the spaxels they fail within the roi
  //______________________________________________________________________
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
    result = Py_BuildValue("(OOOOO)", spaxel_flux_arr, spaxel_weight_arr, spaxel_var_arr,
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
    "cube_match_sky",             /* m_name */
    "find point cloud matches for each spaxel center",  /* m_doc */
    -1,                          /* m_size */
    cube_methods,      /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_cube_match_sky(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
