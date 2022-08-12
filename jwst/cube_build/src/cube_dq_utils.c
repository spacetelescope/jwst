// This is a library of routine to set the DQ plane of the IFU cube

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>

// routines used from cube_utils.c

extern double sh_find_overlap(const double xcenter, const double ycenter, 
			      const double xlength, const double ylength,
			      double xPixelCorner[], double yPixelCorner[]);


//________________________________________________________________________________
// allocate the memory for the spaxel DQ array

int mem_alloc_dq(int nelem, int **idqv) {
    
    const char *msg = "Couldn't allocate memory for output arrays.";

    if (!(*idqv = (int*)calloc(nelem, sizeof(int)))) {
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
			   long ncube, long npt,
			   double *corner1, double *corner2, double *corner3, double *corner4) {
  /* 
     For wavelength plane determine the corners (in xi,eta) of the FOV for MIRI
     Use the 2 extreme slices set by start_region and end_region to define the FOV of the wavelength plane FOV
     Using the min and max coordinates for the extent of the slice on the sky for these two slice - set the
     corners of the FOV. 

     w : wavelength plane
     start_region : starting slice # for channel (slice # in nirspec) 
     end_region : ending slice # for channel     (slice # in nirspec)
     roiw_ave : average roiw for all wavelengths
     zc : array of wavelengths
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

  int slice, c1_use;
  long ipt;
  double wave_distance;
  float c11, c21, c12, c22, length_c1_start, length_c2_start;
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

  // TODO find a better way to find min and max than using hard coded values
  // A little tricky because the min and max is set based on wave_distance & slice #
  float c1_end_min = 10000;
  float c2_end_min = 10000;
  float c1_end_max = -10000;
  float c2_end_max = -10000;

  // Loop over every point cloud member and pull out values
  // 1. Fail withing roiw_ave of wavelength plane
  // and 
  // 2. Are for either of the 2 extreme slices
  
  for (ipt =0; ipt< npt ; ipt++){
    slice = (int)sliceno[ipt];
    wave_distance = fabs(zc[w] - wave[ipt]);

    c11 = -1;  // coord1 for start region
    c21 = -1;  // coord2 for start region
    c12 = -1;  // coord1 for end region
    c22 = -1;  // corrd2 for end region

    if(slice == start_region || slice==end_region){
      // Find all the coordinates on wave slice with slice = start region
      // These points will define corner 1 (min c2) and corner 2 (max c2)
      if (wave_distance < roiw_ave && slice == start_region){ 
	c11 = coord1[ipt];
	c21 = coord2[ipt];
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
      }

      // Find all the coordinates on wave slice with slice = start region
      // These points will define corner 4 (min c2) and corner 3 (max c2)
      if (wave_distance < roiw_ave && slice == end_region){
	c12 = coord1[ipt];
	c22 = coord2[ipt];
	// for the end region find min and max c1,c2
	if (c12 != -1){
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
    }
  } // end looping over point cloud

  // Make sure the 2 extreme slices are found on the FOV. Not finding both can occur for edge wavelength planes
  // or empty wavelength planes between channels

  if (ic1_start_min == -1 || ic1_start_max == -1 || ic1_end_min == -1 || ic1_end_max == -1){
    status = 1;
    return status;
  } else {
  // Find the length in the start and end regions for c1 and c2. This will help define which coordinates to use to set corners.
  // Because we do not know the orientation on the sky pick the  longest length to set how to pick corners. 
    length_c1_start = c1_start_max - c1_start_min;
    length_c2_start = c2_start_max - c2_start_min;

    c1_use = 1; // use the c1 coords to set corners 
    if(length_c1_start < length_c2_start){
      c1_use = 0;   // use the c2 coords to set the corners
    }

    if (c1_use == 0) {
      corner1[0] = coord1[ic2_start_min];
      corner1[1] = coord2[ic2_start_min];
    
      corner2[0] = coord1[ic2_start_max];
      corner2[1] = coord2[ic2_start_max];
    
      corner3[0] = coord1[ic2_end_max];
      corner3[1] = coord2[ic2_end_max];

      corner4[0] = coord1[ic2_end_min];
      corner4[1] = coord2[ic2_end_min];
    } else {
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
// MIRI DQ routine. Find the overlap of the FOV for the wavelength slice in IFU cube

int overlap_fov_with_spaxels(int overlap_partial, int overlap_full,
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
        a. overlap_partial
        b  overlap_full
        bit_wise combination of these values is allowed to account for
        dithered FOVs.

        Parameters
        ----------
	overlap_partial : int
	overlap_full : int
	cdelt1 : double
          IFU cube naxis 1 spatial scale
	cdelt2 : double
          IFU cube naxis 1 spatial scale
	naxis1 : int 
          IFU cube naxis 1 size
	naxis2 : int
          IFU cube naxis 1 size
	xcenters : double array
          IFU naxis 1 array of xi values
	ycenters : double array
          IFU naxis 2 array of eta values 
        xi_corner: xi coordinates of the 4 corners of the FOV on the wavelength plane
        eta_corner: eta coordinates of the 4 corners of the FOV on the wavelength plane

        Sets
        -------
        wave_slice_dq: array containing intermediate dq flag

  */

  // loop over spaxels in the wavelength plane and set slice_dq
  // roughly find the spaxels that might be overlapped

  int i, ixy,ix,iy;
  double area_box, tolerance_dq_overlap, x1, x2, y1, y2, area_overlap, overlap_coverage;
  double ximin = xi_corner[0];
  double etamin = eta_corner[0];
  double ximax = ximin;
  double etamax = etamin;
  for (i=1; i<4 ; i++){
    if (xi_corner[i] < ximin){ ximin = xi_corner[i];}
    if (xi_corner[i] > ximax){ ximax = xi_corner[i];}
    if (eta_corner[i] < etamin){ etamin = eta_corner[i];}
    if (eta_corner[i] < etamax){ etamax = eta_corner[i];}
  }

  area_box = cdelt1 * cdelt2;
  tolerance_dq_overlap = 0.05;  //spaxel has to have 5% overlap to flag in FOV
  // loop over cube xcenters and cube ycenters
  for (ix = 0; ix < naxis1; ix++){
    x1 = (xcenters[ix] - cdelt1)/2;
    x2 = (xcenters[ix] + cdelt1)/2;
    if(x1 > ximin && x2 < ximax){
      for (iy = 0; iy< naxis2; iy++){

	y1 = (ycenters[ix] - cdelt2)/2;
	y2 = (ycenters[ix] + cdelt2)/2;
	if(y1 > etamin && y2 < etamax){
	  ixy = iy * naxis1 + ix;
	  area_overlap = sh_find_overlap(xcenters[ix],
					 ycenters[iy],
					 cdelt1, cdelt2,
					 xi_corner, eta_corner);

	  overlap_coverage = area_overlap / area_box;
	  
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
  
  return 0 ;
}

//________________________________________________________________________________
// Routine to setting NIRSpec dq plane for each wavelength plane

long match_wave_plane_nirspec(double wave_plane,
			      double roiw_ave,
			      double coord1[],
			      double coord2[],
			      double wave[],
			      double sliceno[], 
			      long npt,
			      double *c1_min, double* c2_min,
			      double *c1_max, double *c2_max,
			      int *match_slice){

  /* 
     NIRSpec dq plane is set by mapping each slice to IFU wavelength plane 
     This routine maps each slice to sky and finds the min and max coordinates on the sky
     of the slice. 

     wave_plane : wavelength of current  plane
     slicevalue : slice # 1 to 30 
     roiw_ave : average roiw for all wavelengths
     coord1, coord2: tangent project coordinate of pt cloud
     wave : point cloud wavelength values
     sliceno: slice value of point cloud.
     npt: number of point cloud elements

     return:
     c1_min, c2_min, c1_max, c2_max, match_slice
   */

  long ipt =0;
  double wave_distance;
  double slice;
  long ii = 0;
  int i =0 ; 
  
  // initialize the values
  float minvalue = 10000.0;
  float maxvalue = -10000.0;
  for (i = 0; i < 30; i++){
    c1_min[i] = minvalue;
    c2_min[i] = minvalue;

    c1_max[i] = maxvalue;
    c2_max[i] = maxvalue;

    match_slice[i] = 0;
  }

  for (ipt =0; ipt< npt ; ipt++){
    slice = sliceno[ipt];

    wave_distance = fabs(wave_plane - wave[ipt]);
    double c1 = coord1[ipt];
    double c2 = coord2[ipt];

    // Find all the coordinates that fall on wavelength plane
    if(wave_distance < roiw_ave){

      int islice = (int)slice -1 ;
      
      if (c1< c1_min[islice] ){
	c1_min[islice] = c1;
      }

      if (c2 < c2_min[islice] ){
	c2_min[islice] = c2;
      }
      
      if (c1 > c1_max[islice]){
	c1_max[islice] = c1;
      }
      
      if (c2> c2_max[islice]){
	c2_max[islice] = c2;
      }

      ii = ii + 1 ;
    }
  }
  // find which slices have a c1,c2 min and max found 
  if (ii > 0) {
    for (i = 0; i< 30; i++){
      if (c1_min[i] != minvalue && c2_min[i] != minvalue &&
	  c1_max[i] != maxvalue && c2_max[i] != maxvalue){
	
	if (c1_min[i] != c1_max[i] && c2_min[i] != c2_max[i]){
	  match_slice[i] = 1;
	}
      }
    }
  }
  return ii;

}

//________________________________________________________________________________
// NIRSpec  Find the overlap of the slices  for the wavelength slice with sky
int overlap_slice_with_spaxels(int overlap_partial,
			       double cdelt1, double cdelt2,
			       int naxis1, int naxis2,
			       double xstart, double ystart,
			       double xi_min, double eta_min,
			       double xi_max, double eta_max,
			       int wave_slice_dq[]) {

  /* 
     Set the initial dq plane indicating if the input data falls on a spaxel

     This algorithm assumes the input data falls on a line in the IFU cube, which is
     the case for NIRSpec slices. The NIRSpec slice's endpoints are used to determine
     which IFU spaxels the slice falls on to set an initial dq flag.
     
     Parameters
     ----------
     overlap_partial: intermediate dq flag

     wavelength: the wavelength bin of the IFU cube working with

     Sets
     ----
     wave_slice_dq : array containing intermediate dq flag

     Bresenham's Line Algorithm to find points a line intersects with grid.

     Given the endpoints of a line find the spaxels this line intersects.

     Returns
     -------
     wave_slice_dq filled in 

  */
 
  int error, ystep, y, x, yuse, xuse;
  int index;
  //set up line - convert to integer values
  int x1 = (int)((xi_min - xstart) / cdelt1);
  int y1 = (int)((eta_min - ystart) / cdelt2);
  int x2 = (int)((xi_max - xstart) / cdelt1);
  int y2 = (int)((eta_max - ystart) / cdelt2);

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

  if (x1 > x2){
    x1 = x2;
    x2 = x1;

    y1 = y2;
    y2 = y1;
  }

  // Recalculate differences
  dx = x2 - x1;
  dy = y2 - y1;

  //calculate error
  error = (int)(dx / 2.0);
  ystep = -1;
  if (y1 < y2){
    ystep = 1;
  }

  // iterate over grid to generate points between the start and end of line
  y = y1;

  for (x = x1; x< (x2 + 1); x++){
    yuse  = y;
    xuse = x ;
    if (is_steep){
	yuse = x;
	xuse = y;
      }

    index = (yuse * naxis1) + xuse;

    wave_slice_dq[index] = overlap_partial;
    error -= abs(dy);
    if (error < 0){
      y += ystep;
      error += dx;
    }
  }
  return 0;
}

//________________________________________________________________________________
// Set the spaxel dq = 0. This is used when not determining the FOV on the sky for
// setting the DQ plane. This is case for internalCal type cubes


int set_dqplane_to_zero(int ncube, int **spaxel_dq){

    int *idqv;  // int vector for spaxel
    if (mem_alloc_dq(ncube, &idqv)) return 1;
     *spaxel_dq = idqv;
    return 0;
}

//________________________________________________________________________________
// Main MIRI routine to set DQ plane

int dq_miri(int start_region, int end_region, int overlap_partial, int overlap_full,
	    int nx, int ny, int nz,
	    double cdelt1, double cdelt2, double roiw_ave,
	    double *xc, double *yc, double *zc,
	    double *coord1, double *coord2, double *wave,
	    double *sliceno,
	    long ncube, long npt, 
	    int **spaxel_dq) {

  int status, status_wave, w, nxy, i, istart, iend, in, ii;
  double xi_corner[4], eta_corner[4];
  
  int *idqv ;  // int vector for spaxel
  if (mem_alloc_dq(ncube, &idqv)) return 1;

  double corner1[2];
  double corner2[2];
  double corner3[2];
  double corner4[2];
    
  // for each wavelength plane find the 2 extreme slices to set FOV. Use these two extreme slices to set up the
  // corner of the FOV for each wavelength

  nxy = nx * ny;
  int wave_slice_dq[nxy];
  // Loop over the wavelength planes and set DQ plane 
  for (w = 0; w  < nz; w++) {
    
    for( i = 0; i < nxy; i ++){
      wave_slice_dq[i] = 0;
    }
    status_wave = 0;   
    status_wave =  corner_wave_plane_miri( w, start_region, end_region, roiw_ave, zc,
					  coord1, coord2, wave, sliceno, ncube, npt,
					  corner1, corner2, corner3, corner4);
    if( status_wave == 0){ // found min and max slice on wavelength plane

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
    }
    istart = nxy*w;
    iend = istart + nxy;
    for( in = istart; in < iend; in ++){
      ii = in - istart;
      if(status_wave == 0){
	idqv[in] = wave_slice_dq[ii];
      }else{
	idqv[in] = 0;
      }
    }

  } // end loop over wavelength

  *spaxel_dq = idqv;

  return 0;
}

//________________________________________________________________________________
//Main NIRSpec routine to set up dq plane
 
int dq_nirspec(int overlap_partial,
	       int nx, int ny, int nz,
	       double cdelt1, double cdelt2, double roiw_ave,
	       double *xc, double *yc, double *zc,
	       double *coord1, double *coord2, double *wave,
	       double *sliceno,
	       long ncube, long npt,
	       int **spaxel_dq) {

  /*
    Set an initial DQ flag for the NIRSPEC IFU cube based on FOV of input data.

    Map the FOV of each NIRSpec slice  to the DQ plane and set an initial DQ
    flagging. For  NIRSpec the 30 different slices map to different FOV on the
    range of wavelengths.  The FOV of the slice is really just a line, so instead of using
    the routines the finds the overlap between a polygon and regular grid-
    which is used for MIRI - an algorithm that determines the spaxels that
    the slice line intersects is used instead.
    
    Parameter
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
  
  int w, islice, status, status_wave, nxy, j;
  long istart, in, iend, ii, i;
  double c1_min, c2_min, c1_max, c2_max;
  int *idqv ;  // int vector for spaxel
  idqv = (int*)calloc(ncube, sizeof(int));

  for (i = 0; i< ncube; i++){
    idqv[i] = 0;
  }
  
  nxy = nx * ny;

  for (w = 0; w  < nz; w++) {
    long imatch = 0;
    double c1_min[30];
    double c1_max[30];
    double c2_min[30];
    double c2_max[30];
    int match_slice[30];

    // At each wavelength plane find the min and max of the
    // tangent plane coordinates for each slice
    imatch =  match_wave_plane_nirspec(zc[w], roiw_ave,
				       coord1, coord2, wave,
				       sliceno, npt,
				       c1_min, c2_min,
				       c1_max, c2_max,
				       match_slice);

    int wave_slice_dq[nxy];
    for (j =0; j< nxy; j++){
      wave_slice_dq[j] = 0;
    }

    int slice_found = 0; 

    if( imatch > 0){ // some matches were found on the wavelength slice
      for (islice = 0; islice< 30 ; islice++){
	float slice_c1_min = c1_min[islice];
	float slice_c1_max = c1_max[islice];

	float slice_c2_min = c2_min[islice];
	float slice_c2_max = c2_max[islice];
	
	if( match_slice[islice] == 1){
	  slice_found = 1;

	  // at the wavelength plane find the overlap of each slice on
	  // output spaxel plane

	  float xstart = xc[0];
	  float ystart = yc[0];
	  status = overlap_slice_with_spaxels(overlap_partial,
					      cdelt1, cdelt2,
					      nx, ny,
					      xstart, ystart,
					      slice_c1_min, slice_c2_min,
					      slice_c1_max, slice_c2_max,
					      wave_slice_dq);
	} // end loop if slice has match
      } // end loop over slices 
	
    } // end loop over imatch > 0 match found for wavelength
    if (imatch > 0 && slice_found ==1){
      istart = nxy*w;
      iend = istart + nxy;

      for (in = istart; in < iend; in ++){
	ii = in - istart;
	idqv[in] = wave_slice_dq[ii];
      }
    }
  } // end of wavelength
  *spaxel_dq = idqv;

  return 0;
}

