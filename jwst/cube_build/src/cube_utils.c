
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <stdbool.h>

#define CP_LEFT 0
#define CP_RIGHT 1
#define CP_BOTTOM 2
#define CP_TOP 3


void addpoint (double x, double y, double xnew[], double ynew[], int *nVertices2){
  xnew[*nVertices2] = x;
  ynew[*nVertices2] = y;
  *nVertices2 = *nVertices2 + 1;
}


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


double find_area_poly(int nVertices,double xPixel[],double yPixel[]){
  
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
  double xPixel[9]= {0.0};
  double yPixel[9] = {0.0};
  double xnew[9]= {0.0}; 
  double ynew[9]= {0.0}; 

  // initialize xPixel, yPixel to the detector pixel corners.
  // xPixel,yPixel is become the clipped polygon vertices inside the cube pixel
  for (int i = 0; i < 4; i++) {
    xPixel[i] = xPixelCorner[i];
    yPixel[i] = yPixelCorner[i];
  }
  xPixel[4] = xPixelCorner[0];
  yPixel[4] = yPixelCorner[0];

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
    areaClipped = find_area_poly(nVertices,xPixel,yPixel);
  }
  
  return areaClipped;
}
