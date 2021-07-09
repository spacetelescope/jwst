/*
 
Main function for Python: setup_dq
Python signature: setup_dq(nplane, cdelt1, cdelt2, zcdel3, ....)

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
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_dq_plane_numpy_api    //WHAT IS THIS AND WHERE IS IT USED???
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#define CP_LEFT 0
#define CP_RIGHT 1
#define CP_BOTTOM 2
#define CP_TOP 3

typedef double dbl_type;
static const int npy_dbl = NPY_DOUBLE;


//________________________________________________________________________________
int mem_alloc(long nelem, int **idqv) {
    int  *i;
    const char *msg = "Couldn't allocate memory for output arrays.";

    //printf( " number of elements allocating memory for %lu \n" , nelem);

    // dq
    i = (int*)calloc(nelem, sizeof(int));
    if (i) {
        *idqv = i;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }
    return 0;
}


//________________________________________________________________________________			 
int corner_wave_plane_miri(int w, int start_region, int end_region,
			   dbl_type roiw_ave,
			   dbl_type *zc,
			   dbl_type *coord1, dbl_type *coord2, dbl_type *wave,
			   dbl_type *sliceno,
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
  int ic1_start_min = 0;
  int ic1_start_max = 0;
  int ic2_start_min = 0;
  int ic2_start_max = 0;

  float c1_start_min = 10000;
  float c2_start_min = 10000;
  float c1_start_max = -10000;
  float c2_start_max = -10000;

  int ic1_end_min = 0;
  int ic1_end_max = 0;
  int ic2_end_min = 0;
  int ic2_end_max = 0;

  float c1_end_min = 10000;
  float c2_end_min = 10000;
  float c1_end_max = -10000;
  float c2_end_max = -10000;

  for (int ipt =0; ipt< npt ; ipt++){
    int slice = (int)sliceno[ipt];
    double wave_distance = fabs(zc[w] - wave[ipt]);
	//if(slice == end_region && wave_distance < roiw_ave){
	// printf("wave_distance %f %f %f %f %i %i %i\n", zc[w], wave[ipt],wave_distance, roiw_ave, slice, end_region,ipt);
	//}

    float c11 = -1;  // coord1 for start region
    float c21 = -1;  // coord2 for start region
    float c12 = -1;  // coord1 for end region
    float c22 = -1;  // corrd2 for end region
	
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
  } // end looping over point cloud
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


  return status;
}

//________________________________________________________________________________
int slice_wave_plane_nirspec(int w, int slicevalue,
		      dbl_type roiw_ave,
		      dbl_type *zc,
		      dbl_type *coord1, dbl_type *coord2, dbl_type *wave,
		      dbl_type *sliceno,
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
  double c1 = -1;  // coord1 for slice region
  double c2 = -1;  // coord2 for slice region
 
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
    //printf(" problem finding c1,c2 min and max in wavelength plane, slice # %i %i \n",w, slicevalue);
    status = 1; 
  } else {
    //printf(" found min max of slice for wavelength %f %i %f %f %f %f \n" , zc[w],slicevalue, *c1_min, *c1_max, *c2_min, *c2_max ); 
  }
 
  return status;
}


//_______________________________________________________________________    

int insideWindow(int edge, double x, double y, 
                 double left,double right, double top, double bottom){
        switch(edge)
        {
        case CP_LEFT:
          return (x > left);
        case CP_RIGHT:
          return (x < right);
        case CP_BOTTOM:
          return (y > bottom);
        case CP_TOP:
          return (y < top);
        }
        return 0;
}

int calcCondition(int edge, double x1, double y1, double x2, double y2, 
                  double left, double right, double top, double bottom) {
  int stat1 = insideWindow(edge,x1,y1,left,right,top,bottom);
  int stat2 = insideWindow(edge,x2,y2,left,right,top,bottom);
  if(!stat1 && stat2)   return 1;
  if(stat1  && stat2)   return 2;
  if(stat1  && !stat2)  return 3;
  if(!stat1 && !stat2)  return 4;
  return 0; //never executed

}

void solveIntersection(int edge ,double x1,double y1,double x2,double y2,
                       double *x,double *y,
                       double left, double right, double top, double bottom){
  float m = 0;
  if(x2 != x1) m = ((double)(y2-y1)/(double)(x2-x1));
  switch(edge)
    {
    case CP_LEFT:
      *x = left;
      *y = y1 + m * (*x - x1);
      break;
    case CP_RIGHT:
      *x = right;
      *y = y1 + m * (*x - x1);
      break;
    case CP_BOTTOM:
      *y = bottom;
      if(x1 != x2)
        *x = x1 + (double)(1/m) * (*y - y1);
      else
        *x = x1;
      break;
    case CP_TOP:
      *y = top;
      if(x1 != x2)
        *x = x1 + (double)(1/m) * (*y - y1);
      else
        *x = x1;
      break;
    }
}


void addpoint (double x, double y, double xnew[], double ynew[], int *nVertices2){
  xnew[*nVertices2] = x;
  ynew[*nVertices2] = y;
  *nVertices2++;
}


double findAreaPoly(int nVertices,double xPixel[],double yPixel[]){
  
  double areaPoly = 0.0;
  double xmin = xPixel[0];
  double ymin = yPixel[0];

  for (int i = 1; i < nVertices; i++){ 
    if(xPixel[i] < xmin) xmin = xPixel[i];
    if(yPixel[i] < ymin) ymin = yPixel[i];
  }
  
  for (int i = 0; i < nVertices-1; i++){
    double area = ( xPixel[i]- xmin)*(yPixel[i+1]-ymin) - (xPixel[i+1]-xmin)*(yPixel[i]-ymin);
    areaPoly = areaPoly + area;
  }
  areaPoly = 0.5* areaPoly;
  return fabs(areaPoly);
}
  

//________________________________________________________________________________
double sh_find_overlap(const double xcenter, const double ycenter, 
		      const double xlength, const double ylength,
		      double xPixelCorner[],double yPixelCorner[])
{
  // user the Sutherland_hedgeman Polygon Clipping Algorithm to solve the overlap region
  // first clip the y-z detector plane by the cube's yz rectangle - find the overlap area
  // Then clip the x-y detector plane by the cube's xy rectangle and find the average x lenghth
  //   overlap: overlap vol = area overlap * x lenght overlap


  double areaClipped = 0.0;
  double top = ycenter + 0.5*ylength;
  double bottom = ycenter - 0.5*ylength;
  double left = xcenter - 0.5*xlength;
  double right = xcenter + 0.5*xlength;

  int nVertices = 4; // input detector pixel vertices

  int MaxVertices = 9;
  double xPixel[9];
  double yPixel[9];
  double xnew[9]; 
  double ynew[9]; 

  // initialize xPixel, yPixel to the detector pixel corners.
  // xPixel,yPixel is become the clipped polygon vertices inside the cube pixel
  for (int i = 0; i < 5; i++) {
    xPixel[i] = xPixelCorner[i];
    yPixel[i] = yPixelCorner[i];
  }


  for (int i = 0 ; i < 4; i++) { // 0:left, 1: right, 2: bottom, 3: top
    int nVertices2 = 0;
    for (int j = 0; j< nVertices; j++){
      double x1 = xPixel[j];
      double y1 = yPixel[j];
      double x2 = xPixel[j+1];
      double y2 = yPixel[j+1];

      int condition = calcCondition(i,x1,y1,x2,y2,
                                    left,right,top,bottom);

      double x = 0;
      double y = 0;

      switch(condition)
        {
        case 1:
          solveIntersection(i,x1,y1,x2,y2,
                            &x, &y,
                            left,right,top,bottom);
          
          
          addpoint (x, y, xnew, ynew, &nVertices2);
          addpoint (x2, y2, xnew, ynew, &nVertices2);
          break;
        case 2:
          addpoint (x2, y2, xnew, ynew, &nVertices2);
          break;
        case 3:
          solveIntersection(i,x1,y1,x2,y2,
                            &x, &y,
                            left,right,top,bottom);
          addpoint (x, y, xnew, ynew, &nVertices2);
          break;
        case 4:
          break;
        }
    }// loop over j  corners


   addpoint (xnew[0], ynew[0], xnew, ynew, &nVertices2); // closed polygon
    

    if(nVertices2 >  MaxVertices ) {
      printf( " sh_find_overlap:: failure in finding the clipped polygon, nVertices2 > 9 \n ");
      exit(EXIT_FAILURE);
    }
    nVertices = nVertices2-1;

    for (int k = 0; k< nVertices2; k++){
      xPixel[k] = xnew[k];
      yPixel[k] = ynew[k];
    }

    //update 

  } // loop over top,bottom,left,right
  nVertices++;
  if(nVertices > 0) {
    areaClipped = findAreaPoly(nVertices,xPixel,yPixel);
  }
  
  return areaClipped;
}


//________________________________________________________________________________
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
  

   // # set for a single wavelength
//else:
//spaxel_dq[wmin, :] = np.bitwise_or(spaxel_dq[wmin, :],
//                                           wave_slice_dq)
					   
}


//________________________________________________________________________________
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
	    dbl_type cdelt1, dbl_type cdelt2, dbl_type roiw_ave,
    dbl_type *xc, dbl_type *yc, dbl_type *zc,
	    dbl_type *coord1, dbl_type *coord2, dbl_type *wave,
	    dbl_type *sliceno,
	    long ncube, int npt, int imin, int imax,
	    int **spaxel_dq) {

  int *idqv = NULL;  // int vector for spaxel

    if (mem_alloc(ncube, &idqv)) return 1;

    // Set all data to zero
    printf( " number of elements allocating memory for %lu \n" , ncube);
    for (long i = 0; i < ncube; i++){
      idqv[i] = 0;
    }

    double corner1[2];
    double corner2[2];
    double corner3[2];
    double corner4[2];
    
    //printf(" imin imax %i %i \n ", imin, imax);
    // for each wavelength plane find the 2 extreme slices to set FOV. Use these two extreme slices to set up the
    // corner of the FOV for each wavelength

    int nxy = nx * ny;
    
    for (int w = imin; w  < imax; w++) {

      int status =  corner_wave_plane_miri( w, start_region, end_region, roiw_ave, zc,
				   coord1, coord2, wave, sliceno, ncube, npt,
				   corner1, corner2, corner3, corner4);

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
	

      
      //printf(" corner %f %f %f %f %f %f %f %f for wavelength %i  \n", corner1[0], corner1[1],corner2[0], corner2[1],
      //	     corner3[0], corner3[1],corner4[0], corner4[1],w);

      int wave_slice_dq[nxy];
      for( int i = 0; i < nxy; i ++){
	wave_slice_dq[i] = 0;
      }
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
	  if(in == 0){
	    printf(" found match %i \n",idqv[in]);
	  }
	}
    } // end loop over wavelength

    *spaxel_dq = idqv;
    free(idqv);

    return 0;
}


//________________________________________________________________________________
// set up the dq plane for NIRSPEC
 
int dq_nirspec(int overlap_partial,
	    int nx, int ny, int nz,
	    dbl_type cdelt1, dbl_type cdelt2, dbl_type roiw_ave,
	    dbl_type *xc, dbl_type *yc, dbl_type *zc,
	    dbl_type *coord1, dbl_type *coord2, dbl_type *wave,
	    dbl_type *sliceno,
	    long ncube, int npt, int imin, int imax,
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

  if (mem_alloc(ncube, &idqv)) return 1;

  // Set all data to zero
  //printf( " number of elements allocating memory for %lu \n" , ncube);
  for (long i = 0; i < ncube; i++){
    idqv[i] = 0;
  }

  //  for each of the 30 slices - find the projection of this slice
  //     onto each of the IFU wavelength planes.


  for (int w = imin; w  < imax; w++) {
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
    if (n_slice_found ==0){
      printf(" no slices found on wavelength plane %i %f \n",w, zc[w]);
    }
	
  } // end of wavelength
  *spaxel_dq = idqv;
  free(idqv);

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


static PyObject * setup_dq(PyObject *module, PyObject *args) {
  PyObject *result = NULL, *xco, *yco, *zco, *coord1o, *coord2o, *waveo, *slicenoo;
  
  double cdelt1, cdelt2, roiw_ave; 
  int start_region, end_region, overlap_partial, overlap_full, nx, ny, nz, imin, imax, npt;
  long ncube;
  int instrument;
  int *spaxel_dq=NULL;
  int free_xc=0, free_yc=0, free_zc=0, free_coord1=0, free_coord2 =0, free_wave=0, free_sliceno=0, status=0;
  
  PyArrayObject *xc, *yc, *zc, *coord1, *coord2, *wave, *sliceno;
  PyArrayObject *spaxel_dq_arr=NULL;
  npy_intp npy_ncube = 0;


  if (!PyArg_ParseTuple(args, "iiiiiiiidddOOOOOOOliii:setup_dq",
			&instrument,
			&start_region, &end_region, &overlap_partial, &overlap_full,
			&nx, & ny, &nz,
			&cdelt1, &cdelt2, &roiw_ave,
			&xco, &yco, &zco,
			&coord1o, &coord2o, &waveo,&slicenoo,
			&ncube, &npt, &imin, &imax)) {
    return NULL;
  }


  // check that input parameters are valid:
  if (nx < 0) {
    PyErr_SetString(PyExc_ValueError,
		    "'nx' must be a strictly positive number.");
    return NULL;
  }

  if (ny < 0) {
    PyErr_SetString(PyExc_ValueError,
		    "'ny' must be a strictly positive number.");
    return NULL;
  }


  if ((cdelt1 < 0) || (cdelt2 < 0)) {
    PyErr_SetString(PyExc_ValueError,
		    "'cdelt1' and 'cdelt2' must be a strictly positive number.");
    return NULL;
  }

  if ((imin < 0) || (imax < 0)) {
    PyErr_SetString(PyExc_ValueError,
		    "'imin' and 'imax' must be a strictly positive number.");
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
      (!(sliceno = ensure_array(slicenoo, &free_sliceno))))

    {
      goto cleanup;
    }


  int ncx = (int) PyArray_Size((PyObject *) coord1);
  int ncy = (int) PyArray_Size((PyObject *) coord2);
  if (ncx != ncy) {
    PyErr_SetString(PyExc_ValueError,
		    "Input coordinate arrays of unequal size.");
    goto cleanup;
  }
  
  int ncz = (int) PyArray_Size((PyObject *) wave);
  if (!ncx  || !ncz) {
    // 0-length input arrays. Nothing to clip. Return 0-length arrays

    spaxel_dq = (PyArrayObject*) PyArray_EMPTY(1, &npy_ncube, NPY_INT, 0);
    if (!spaxel_dq_arr) goto fail;


    result = Py_BuildValue("(O)", spaxel_dq_arr);
	
    goto cleanup;
  }

  if (instrument == 0){
    status = dq_miri(start_region, end_region,overlap_partial, overlap_full,
		     nx,ny,nz,
		     cdelt1, cdelt2, roiw_ave,
		     (dbl_type *) PyArray_DATA(xc),
		     (dbl_type *) PyArray_DATA(yc),
		     (dbl_type *) PyArray_DATA(zc),
		     (dbl_type *) PyArray_DATA(coord1),
		     (dbl_type *) PyArray_DATA(coord2),
		     (dbl_type *) PyArray_DATA(wave),
		     (dbl_type *) PyArray_DATA(sliceno),
		     ncube, npt, imin, imax,
		     &spaxel_dq);

  } else{
    status = dq_nirspec(overlap_partial,
			nx,ny,nz,
			cdelt1, cdelt2, roiw_ave,
			(dbl_type *) PyArray_DATA(xc),
			(dbl_type *) PyArray_DATA(yc),
			(dbl_type *) PyArray_DATA(zc),
			(dbl_type *) PyArray_DATA(coord1),
			(dbl_type *) PyArray_DATA(coord2),
			(dbl_type *) PyArray_DATA(wave),
			(dbl_type *) PyArray_DATA(sliceno),
			ncube, npt, imin, imax,
			&spaxel_dq);


  }
  if (status) {
    goto fail;

  } else {
    // create return tuple:
    npy_ncube = (npy_ulong) ncube;
    
    spaxel_dq_arr = (PyArrayObject*) PyArray_SimpleNewFromData(1, &npy_ncube, NPY_INT, spaxel_dq);
								 
    if (!spaxel_dq_arr) goto fail;
    spaxel_dq = NULL;


    PyArray_ENABLEFLAGS(spaxel_dq_arr, NPY_ARRAY_OWNDATA);
    
    result = Py_BuildValue("(O)", spaxel_dq_arr);
    goto cleanup;
  }

 fail:
  Py_XDECREF(spaxel_dq_arr);
  free(spaxel_dq);
    

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
  if (free_sliceno) Py_XDECREF(sliceno);
  return result;
}


static PyMethodDef dq_methods[] =
{
    {
        "setup_dq",  setup_dq, METH_VARARGS,
        "setup_dq(put in doc string)"
    },
    {0, 0}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "dq_plane",                  /* m_name */
    "set up dq plane",           /* m_doc */
    -1,                          /* m_size */
    dq_methods,                  /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_dq_plane(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
