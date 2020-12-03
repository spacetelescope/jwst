/* 
NAME:

  POLYCLIP 

DESCRIPTION:

  Recursively clips an input polygon or set of polygons to a range of
  pixels on a square grid using Sutherland-Hodgeman clipping,
  returning the clipped polygons (for a single input polygon only),
  along with the clipped polygon areas.  Uses REVERSE_INDICES-like
  vectors (see HISTOGRAM) to permit decoding the return vectors into
  individual polygons clipped for each pixel.  Much faster (~50x) than
  the IDL-loop version polyclip.pro.  

REFERENCE: 

Sutherland, I, and Hodgman, G, "Reentrant Polygon Clipping", Graphics
  and Image Processing, Vol 17, p32, Jan, 1974

CALLING SEQUENCE:

         tmp=call_external(polyclip_path,'polyclip_single',$
                          VALUE= $
                          [1b, 1b,   1b,    1b,  $
                           0b,0b,1b, $
                           0b, $
                           0b, $
                           0b, $
                           0b,0b,0b], $
                          left,right,bottom,top, $ ; bound-box to consider
                          px,py,nverts, $          ; input polygon indices
                          inds, $                  ; OUT: x,y inds within array
                          nclip_poly, $       ; OUT: final number clipped polys
                          areas, $            ; OUT: output areas
                          px_out,py_out,ri_out) ; OUT: clipped poly verts
        
     or, for multiple polygons at once:

        tmp=call_external(polyclip_path,'polyclip_multi',$
                          VALUE= $
                          [0b,0b,0b,0b,$
			   0b,0b, $
                           1b,   0b, $
                           0b, $
                           0b, $
                           0b], $
                          left,right,bottom,top, $ ; bounding arrays (n_poly)
                          px,py, $            ; input polygon indices
                          n_poly,poly_inds, $ ; indices per poly into px,py
                          inds, $             ; OUT: x,y inds within array
                          nclip_poly, $       ; OUT: final number clipped polys
                          areas)              ; OUT: output areas

INPUTS: 

  left,right,bottom,top: The bounding-box range of pixel coordinates
    to clip against.  For multiple input polygons, these bounding
    boxes are arrays of length n_poly.  The total number of pixels
    enclosed by the bounding box(es) is npix.

  px, py: The polygon(s) as a series of points (not closed).

  == POLYCLIP_SINGLE only:

  nverts: The number of vertices in the polygon.

  == POLYCLIP_MULTI only:

  n_poly: The number of polygons passed in.

  poly_inds: A reverse index style vector, as long integer array
    pre-initialized to a length of n_poly+1.
    poly_inds[i]:poly_inds[i+1]-1 gives the range of indices within px
    and py corresponding to input polygon i [0..n_poly-1].  On output,
    see below.


OUTPUTS:

  inds: a 2xn array giving the X,Y coordinates of pixels corresponding
    to output polygons and areas.  Indexed by poly_inds for the
    multiple case.
  
  nclip_poly: The total number of resulting clipped polygons.

  areas: The output areas of each polygon, must be pre-initialized as
    a vector of length npix.  Indexed by poly_inds for the multiple
    case.

  == POLYCLIP_SINGLE only:

  px_out,py_out: The output polygons, concatenated together.  Must be
    pre-initialized to contain at least npix*(nverts+4) points (for
    convex polygons), or more for general polygons.

  ri_out: The reverse index vector for px_out/py_out, which must be
    pre-initialzed as a long integer vector of length npix+1,
    subsequent entries give the range of indices in the px_out, py_out
    vectors corresponding to that pixel's clipped polygon coordinates.

  == POLYCLIP_MULTIPLE only:
  
  poly_inds: On output, in contains the reverse index vector for inds
    and areas, such that poly_inds[i]:poly_inds[i+1]-1 contains the
    indices into these vectors corresponding to input polygon i.  The
    length will be nclip_poly+1.  Note that the clipped polygons
    themselves cannot be output by POLYCLIP_MULTIPLE.  Also serves as
    an input argument (see above).


EXAMPLE:

  make_dll,'polyclip','polyclip',['polyclip_single','polyclip_multi']
  xy=[76.65864,78.83240,44.87576,79.87564,$
      44.84472,78.90629,76.62761,77.86305]
  px=xy[indgen(4)*2] & py=xy[indgen(4)*2+1] & sz=[128,128]
  plot,[px,px[0]],[py,py[0]],yrange=[77,80]
  for r=!X.CRANGE[0],!X.CRANGE[1],1 do plots,r,!Y.CRANGE
  for r=!Y.CRANGE[0],!Y.CRANGE[1],1 do plots,!X.CRANGE,r
  left=floor(min(px,max=maxx))>0 &   right=floor(maxx)<(sz[0]-1)
  bottom=floor(min(py,max=maxy))>0 & top=floor(maxy)<(sz[1]-1)
  nx=right-left+1 & ny=top-bottom+1
  inds=lonarr(2,nx*ny)
  areas=fltarr(nx*ny*(n_elements(px)+4))
  ri_out=lonarr(nx*ny+1)
  nclip_poly=0L
  tmp=call_external(polyclip_path,'polyclip_single', $
                    VALUE= $
		    [1b, 1b,   1b,    1b,  $
		    0b,0b,1b, $
		    0b, $
		    0b, $
		    0b, $
		    0b,0b,0b], $
		    left,right,bottom,top, $
		    px,py,nverts, $         
		    inds, $                  ; OUT: x,y inds within array
		    nclip_poly, $       ; OUT: final number clipped polys
		    areas, $            ; OUT: output areas
		    px_out,py_out,ri_out) x_out
		    
#############################################################################
 LICENSE

  Copyright (C) 2001,2002,2003,2006, 2007 J.D. Smith

  This file is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published
  by the Free Software Foundation; either version 2, or (at your
  option) any later version.
  
  This file is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this file; see the file COPYING.  If not, write to the
  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA 02110-1301, USA.
#############################################################################*/
/*    $Id$*/

#include <stdlib.h>
#include <stdio.h>

int  polyclip(float *,float *, int, int, int, float *, float *);
void  polyclip_shclip(float, float, int, int, int);
void polyclip_shclose(int, int, int);
int  polyclip_inside(float, float, int, int, int);
void polyclip_intersect(float, float, int, int, int);
float polyclip_area(float *px,float *py, int n);
char polyclip_test();
//void polyclip_test(int , int , int b, int t, float *px, float *py, int nverts, int *inds, int *nclip_poly, float *areas,
//  float *px_out, float *py_out, int *ri_out);

/* Return a unique version number to test compilation success */
char polyclip_test() {
  return 44;
}

/* Clip a single polygon, with polygon output */
void polyclip_single(int argc,void * argv[]) {
  /* polyclip_single(l,r,b,t,px,py,nverts,inds,nclip_poly,areas, $ */
  /*                 px_out,py_out,ri_out) */
  int i,j,l,r,b,t,nverts,nv_clip,indx;
  float *px,*py,*px_out,*py_out,*areas,area;
  int *inds,*nclip_poly,*ri_out;
  
  /* Input */
  l=(int)argv[0]; r=(int)argv[1]; b=(int)argv[2]; t=(int)argv[3];
  px=(float *)argv[4]; py=(float *)argv[5]; nverts=(int)argv[6];

  /* Output */
  inds=(int *)argv[7]; nclip_poly=(int *)argv[8];
  areas=(float *)argv[9];
  px_out=(float *)argv[10]; py_out=(float *)argv[11]; ri_out=(int *)argv[12];

  ri_out[0]=0;
  for(indx=0,i=l;i<=r;i++) {
    for(j=b;j<=t;j++) {
      if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
	area=polyclip_area(px_out,py_out,nv_clip);
	if (area==0.0) continue;
	areas[indx]=area;	/* Discard degenerates */
	(*nclip_poly)++;
	ri_out[indx+1]=ri_out[indx]+nv_clip;
	px_out+=nv_clip; py_out+=nv_clip; /* Offset for next output poly */
	inds[2*indx]=i; inds[2*indx+1]=j;
	indx++;
      } 
    }
  }
}

/* Clip a single polygon, with polygon output */
void polyclip_single2(int l, int r, int b, int t, float *px, float *py, int nverts, int *inds, int *nclip_poly, float *areas,
  float *px_out, float *py_out, int *ri_out) {
  /* polyclip_single(l,r,b,t,px,py,nverts,inds,nclip_poly,areas, $ */
  /*                 px_out,py_out,ri_out) */
  //int i,j,l,r,b,t,nverts
  int i,j,nv_clip,indx;
  float area;
  //float *px,*py,*px_out,*py_out,*areas,area;
  //int *inds,*nclip_poly,*ri_out;
  
  /* Input */
  //l=(int)argv[0]; r=(int)argv[1]; b=(int)argv[2]; t=(int)argv[3];
  //px=(float *)argv[4]; py=(float *)argv[5]; nverts=(int)argv[6];

  /* Output */
  //inds=(int *)argv[7]; nclip_poly=(int *)argv[8];
  //areas=(float *)argv[9];
  //px_out=(float *)argv[10]; py_out=(float *)argv[11]; ri_out=(int *)argv[12];

  ri_out[0]=0;
  for(indx=0,i=l;i<=r;i++) {
    for(j=b;j<=t;j++) {
      if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
  area=polyclip_area(px_out,py_out,nv_clip);
  if (area==0.0) continue;
  areas[indx]=area; /* Discard degenerates */
  
  (*nclip_poly)++;
  ri_out[indx+1]=ri_out[indx]+nv_clip;
  px_out+=nv_clip; py_out+=nv_clip; /* Offset for next output poly */
  inds[2*indx]=i; inds[2*indx+1]=j;
  indx++;
      } 
    }
  }
}
/* Clip multiple polygons (without any output polygons) */
void polyclip_multi(int argc, void* argv[]) {
  /* polyclip_multi(px,py,n_poly,poly_inds, $ */
  /*                inds,nclip_poly,areas)            */
  
  int i,j,k,nv_clip,indx,*l,*r,*b,*t;
  float *px,*py,*px_out,*py_out,*areas,area;
  int n_poly;
  unsigned int *poly_inds;
  int *inds,*nclip_poly, nverts, this_nclip_poly, prev_pind, nv_max;

  /* Input */
  l=(int *)argv[0]; r=(int *)argv[1]; b=(int *)argv[2]; t=(int *)argv[3];
  px=(float *)argv[4]; py=(float *)argv[5]; 
  n_poly=(int)argv[6];
  poly_inds=(unsigned int *)argv[7]; /* poly_inds Input/Output */

  /* Output */
  inds=(int *)argv[8]; nclip_poly=(int *)argv[9];
  areas=(float *)argv[10];
  
  /* Maximal output polygon: input + 4 vertices */
  for(nv_max=0, k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-poly_inds[k];
    if(nverts>nv_max) nv_max=nverts;
  }
  nv_max+=24;			/* for a margin of safety, we include
				   24 more */
  px_out=(float *)malloc((nv_max)*sizeof(float));
  py_out=(float *)malloc((nv_max)*sizeof(float));
  /* Clip each polygon and accumulate results */
  for(indx=0,prev_pind=0,k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-prev_pind;
    this_nclip_poly=0;
    for(i=l[k];i<=r[k];i++) {
      for(j=b[k];j<=t[k];j++) {
	if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
	  area=polyclip_area(px_out,py_out,nv_clip);
	  if (area==0.0) continue; /* Discard degenerates */
	  areas[indx]=area;	
	  this_nclip_poly++;
	  inds[2*indx]=i; inds[2*indx+1]=j;
	  indx++;
	}
      }
    }
    (*nclip_poly)+=this_nclip_poly; /* Number of resulting polygons */
    prev_pind=poly_inds[k+1]; /* Reusing poly_inds as input and output */
    poly_inds[k+1]=poly_inds[k]+this_nclip_poly; /* Reverse index */
    px+=nverts; py+=nverts;	/* Offset to next input poly */
  }
  free(px_out); free(py_out);
}

void polyclip_multi2(int *l, int *r, int *b, int *t, float *px, float *py, int n_poly, int *poly_inds, int *inds, int *nclip_poly, float *areas) {
  /* polyclip_multi(px,py,n_poly,poly_inds, $ */
  /*                inds,nclip_poly,areas)            */
  
  int i,j,k,nv_clip,indx;
  float *px_out,*py_out,area;
  int nverts, this_nclip_poly, prev_pind, nv_max;


  //int i,j,k,nv_clip,indx,*l,*r,*b,*t;
  //float *px,*py,*px_out,*py_out,*areas,area;
  //int n_poly;
  //unsigned int *poly_inds;
  //int *inds,*nclip_poly, nverts, this_nclip_poly, prev_pind, nv_max;

  /* Input */
  //l=(int *)argv[0]; r=(int *)argv[1]; b=(int *)argv[2]; t=(int *)argv[3];
  //px=(float *)argv[4]; py=(float *)argv[5]; 
  //n_poly=(int)argv[6];
  //poly_inds=(unsigned int *)argv[7]; /* poly_inds Input/Output */

  /* Output */
  //inds=(int *)argv[8]; nclip_poly=(int *)argv[9];
  //areas=(float *)argv[10];
  //printf("n_poly:%d\n",n_poly);fflush(stdout);
  /* Maximal output polygon: input + 4 vertices */
  for(nv_max=0, k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-poly_inds[k];
      //printf("nverts:%d\n",nverts);fflush(stdout);

    if(nverts>nv_max) nv_max=nverts;
  }
  nv_max+=24;     /* for a margin of safety, we include
           24 more */
  px_out=(float *)malloc((nv_max)*sizeof(float));
  py_out=(float *)malloc((nv_max)*sizeof(float));
  /* Clip each polygon and accumulate results */
  for(indx=0,prev_pind=0,k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-prev_pind;
      //printf("nverts 2:%d\n",nverts);fflush(stdout);

    this_nclip_poly=0;
    for(i=l[k];i<=r[k];i++) {
      for(j=b[k];j<=t[k];j++) {

        //printf("%d %d %f %f %f %f\n",i,j,px[0],px[1],px[2],px[3]);fflush(stdout);
        //printf("%d %d %f %f %f %f\n",i,j,py[0],py[1],py[2],py[3]);fflush(stdout);

  if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
    area=polyclip_area(px_out,py_out,nv_clip);
    //printf("area:%f\n",area);fflush(stdout);
    if (area==0.0) continue; /* Discard degenerates */
    areas[indx]=area; 
    this_nclip_poly++;
    inds[2*indx]=i; inds[2*indx+1]=j;
    indx++;
  }
      }
    }
    (*nclip_poly)+=this_nclip_poly; /* Number of resulting polygons */
    prev_pind=poly_inds[k+1]; /* Reusing poly_inds as input and output */
    poly_inds[k+1]=poly_inds[k]+this_nclip_poly; /* Reverse index */
    px+=nverts; py+=nverts; /* Offset to next input poly */
  }
  free(px_out); free(py_out);
}

void polyclip_multi3(int *l, int *r, int *b, int *t, float *px, float *py, int n_poly, int *poly_inds, int *inds, int *nclip_poly, float *areas, int *index) {
  /* polyclip_multi(px,py,n_poly,poly_inds, $ */
  /*                inds,nclip_poly,areas)            */
  
  int i,j,k,nv_clip,indx;
  float *px_out,*py_out,area;
  int nverts, this_nclip_poly, prev_pind, nv_max;


  //int i,j,k,nv_clip,indx,*l,*r,*b,*t;
  //float *px,*py,*px_out,*py_out,*areas,area;
  //int n_poly;
  //unsigned int *poly_inds;
  //int *inds,*nclip_poly, nverts, this_nclip_poly, prev_pind, nv_max;

  /* Input */
  //l=(int *)argv[0]; r=(int *)argv[1]; b=(int *)argv[2]; t=(int *)argv[3];
  //px=(float *)argv[4]; py=(float *)argv[5]; 
  //n_poly=(int)argv[6];
  //poly_inds=(unsigned int *)argv[7]; /* poly_inds Input/Output */

  /* Output */
  //inds=(int *)argv[8]; nclip_poly=(int *)argv[9];
  //areas=(float *)argv[10];
  //printf("n_poly:%d\n",n_poly);fflush(stdout);
  /* Maximal output polygon: input + 4 vertices */
  for(nv_max=0, k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-poly_inds[k];
      //printf("nverts:%d\n",nverts);fflush(stdout);

    if(nverts>nv_max) nv_max=nverts;
  }
  nv_max+=24;     /* for a margin of safety, we include
           24 more */
  px_out=(float *)malloc((nv_max)*sizeof(float));
  py_out=(float *)malloc((nv_max)*sizeof(float));
  /* Clip each polygon and accumulate results */
  for(indx=0,prev_pind=0,k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-prev_pind;
      //printf("nverts 2:%d\n",nverts);fflush(stdout);

    this_nclip_poly=0;
    for(i=l[k];i<=r[k];i++) {
      for(j=b[k];j<=t[k];j++) {

        //printf("%d %d %f %f %f %f\n",i,j,px[0],px[1],px[2],px[3]);fflush(stdout);
        //printf("%d %d %f %f %f %f\n",i,j,py[0],py[1],py[2],py[3]);fflush(stdout);

  if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
    area=polyclip_area(px_out,py_out,nv_clip);
    //printf("area:%f\n",area);fflush(stdout);
    if (area==0.0) continue; /* Discard degenerates */
    areas[indx]=area; 
    this_nclip_poly++;
    inds[2*indx]=i; inds[2*indx+1]=j;
    index[indx] = k;
    indx++;
  }
      }
    }
    (*nclip_poly)+=this_nclip_poly; /* Number of resulting polygons */
    prev_pind=poly_inds[k+1]; /* Reusing poly_inds as input and output */
    poly_inds[k+1]=poly_inds[k]+this_nclip_poly; /* Reverse index */
    px+=nverts; py+=nverts; /* Offset to next input poly */
  }
  free(px_out); free(py_out);
}

void polyclip_multi4(int *l, int *r, int *b, int *t, float *px, float *py, int n_poly, int *poly_inds, int *x, int*y, int *nclip_poly, float *areas, int *index) {
  /* polyclip_multi(px,py,n_poly,poly_inds, $ */
  /*                inds,nclip_poly,areas)            */
  
  int i,j,k,nv_clip,indx;
  float *px_out,*py_out,area;
  int nverts, this_nclip_poly, prev_pind, nv_max;


  //int i,j,k,nv_clip,indx,*l,*r,*b,*t;
  //float *px,*py,*px_out,*py_out,*areas,area;
  //int n_poly;
  //unsigned int *poly_inds;
  //int *inds,*nclip_poly, nverts, this_nclip_poly, prev_pind, nv_max;

  /* Input */
  //l=(int *)argv[0]; r=(int *)argv[1]; b=(int *)argv[2]; t=(int *)argv[3];
  //px=(float *)argv[4]; py=(float *)argv[5]; 
  //n_poly=(int)argv[6];
  //poly_inds=(unsigned int *)argv[7]; /* poly_inds Input/Output */

  /* Output */
  //inds=(int *)argv[8]; nclip_poly=(int *)argv[9];
  //areas=(float *)argv[10];
  //printf("n_poly:%d\n",n_poly);fflush(stdout);
  /* Maximal output polygon: input + 4 vertices */
  for(nv_max=0, k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-poly_inds[k];
      //printf("nverts:%d\n",nverts);fflush(stdout);

    if(nverts>nv_max) nv_max=nverts;
  }
  nv_max+=24;     /* for a margin of safety, we include
           24 more */
  px_out=(float *)malloc((nv_max)*sizeof(float));
  py_out=(float *)malloc((nv_max)*sizeof(float));
  /* Clip each polygon and accumulate results */
  for(indx=0,prev_pind=0,k=0;k<n_poly;k++) {
    nverts=poly_inds[k+1]-prev_pind;
      //printf("nverts 2:%d\n",nverts);fflush(stdout);

    this_nclip_poly=0;
    for(i=l[k];i<=r[k];i++) {
      for(j=b[k];j<=t[k];j++) {

        //printf("%d %d %f %f %f %f\n",i,j,px[0],px[1],px[2],px[3]);fflush(stdout);
        //printf("%d %d %f %f %f %f\n",i,j,py[0],py[1],py[2],py[3]);fflush(stdout);

  if((nv_clip=polyclip(px,py,nverts,i,j,px_out,py_out))) {
    area=polyclip_area(px_out,py_out,nv_clip);
    //printf("area:%f\n",area);fflush(stdout);
    if (area==0.0) continue; /* Discard degenerates */
    areas[indx]=area; 
    this_nclip_poly++;
    //inds[2*indx]=i; inds[2*indx+1]=j;
    x[indx] = i;
    y[indx] = j;
    index[indx] = k;
    indx++;
  }
      }
    }
    (*nclip_poly)+=this_nclip_poly; /* Number of resulting polygons */
    prev_pind=poly_inds[k+1]; /* Reusing poly_inds as input and output */
    poly_inds[k+1]=poly_inds[k]+this_nclip_poly; /* Reverse index */
    px+=nverts; py+=nverts; /* Offset to next input poly */
  }
  free(px_out); free(py_out);
}
//------------------------------------------------------------------------
// Sutherland-Hodgman clipper code
//------------------------------------------------------------------------

#define LEFT 0
#define RIGHT 1
#define TOP 2
#define BOTTOM 3
#define DONE 4

int in_last[4], first[4]; /* Flags for first and inside, for each side */
float *px_clip,*py_clip; /* pointers for depositing output vertices */
float F[4][2],S[4][2],I[2];	/* Last point X, Y in poly */
int pind;			/* Counter for accumulating output */

int polyclip(float *px, float *py, int n, int i, int j, 
	     float *px_out, float *py_out) {
  int l;
  pind=0; px_clip=px_out; py_clip=py_out;

#ifdef DEBUG
  for(l=0;l<n;l++) printf("%8.5f %8.5f\n",px[l],py[l]);
#endif

  for(l=0;l<4;l++) first[l]=1;
  for(l=0;l<n;l++) 
    polyclip_shclip(px[l],py[l],i,j,LEFT);
  polyclip_shclose(i,j,LEFT);	/* close first->last */
  return pind;
}

/* Reentrant Sutherland-Hodgman Clipper */
/* Recursively clip a polygon with all 4 boundaries of pixel (i,j) */
void polyclip_shclip(float px, float py, int i, int j, int side) {
  int in_p;

#ifdef DEBUG
  if (side < DONE) 
    printf("=== Clipping (%4.2f,%4.2f) pixel %d %d %s\n",px,py,i,j,
	   (side==LEFT)?"LEFT":((side==RIGHT)?"RIGHT":
				((side==TOP)?"TOP":"BOTTOM")));
#endif

  if (side==DONE) { 			/* Done, store the point */
    px_clip[pind]=px; py_clip[pind++]=py;
#ifdef DEBUG
    printf("Added: %f %f\n",px,py);
#endif

    return;
  }
  
  in_p=polyclip_inside(px,py,i,j,side);

  if(first[side]) {
    first[side]=0;
    F[side][0]=px; F[side][1]=py; /* P -> F */
  } else if(in_last[side]^in_p) {	/* Crossed -- compute intersection */
    polyclip_intersect(px,py,i,j,side);
#ifdef DEBUG
    printf("Intersec (%4.2f,%4.2f) -> (%4.2f,%4.2f) => (%4.2f,%4.2f) %s-%s\n",
	   S[side][0],S[side][1],px,py,I[0],I[1],in_last[side]?"in":"out",
	   in_p?"in":"out");
#endif
    polyclip_shclip(I[0],I[1],i,j,side+1); /* Pass this point to the next */
  }
    
  S[side][0]=px; S[side][1]=py;  /* P -> S */
  in_last[side]=in_p;		 /* Save last inside flag */
  if(in_p) polyclip_shclip(px,py,i,j,side+1); 
}

void polyclip_shclose(int i, int j, int side) {
#ifdef DEBUG
  if(side<DONE) 
    printf("Closing pixel %d %d (inlast: %d, F: %7.4f, %7.4f, first: %d) %s\n",
	   i,j,in_last[side],F[side][0],F[side][1],first[side],
	   (side==LEFT)?"LEFT":((side==RIGHT)?"RIGHT":
				((side==TOP)?"TOP":"BOTTOM")));
#endif
  if (side<DONE) {
    if(!first[side]) {
      if(in_last[side]^polyclip_inside(F[side][0],F[side][1],i,j,side)) {
	polyclip_intersect(F[side][0],F[side][1],i,j,side);
	
#ifdef DEBUG
	printf("Intersec (%4.2f,%4.2f) -> (%4.2f,%4.2f) => (%4.2f,%4.2f) last %s\n",
	       S[side][0],S[side][1],F[side][0],F[side][1],I[0],I[1],in_last[side]?"in":"out");
#endif
	
	polyclip_shclip(I[0],I[1],i,j,side+1);
      }
      first[side]=1;
    }
    polyclip_shclose(i,j,side+1);
  }
}
      
int polyclip_inside(float px, float py, int i, int j, int side) {
  switch(side) { 		/* See if inside the edge */
  case LEFT: 
    return (px>=i);
  case RIGHT: 
    return (px<=i+1);
  case TOP: 
    return (py<=j+1);
  case BOTTOM: 
    return (py>=j);
  }
  return -1;
}

void polyclip_intersect(float px, float py,int i, int j, int side) { 
  switch(side) {
  case LEFT:
    I[0]=i;
    I[1]=S[side][1]+(py-S[side][1])/(px-S[side][0])*(i-S[side][0]);
    break;
  case RIGHT: 
    I[0]=i+1;
    I[1]=S[side][1]+(py-S[side][1])/(px-S[side][0])*(i+1-S[side][0]);
    break;
  case TOP: 
    I[0]=S[side][0]+(px-S[side][0])/(py-S[side][1])*(j+1-S[side][1]);
    I[1]=j+1;
    break;
  case BOTTOM:
    I[0]=S[side][0]+(px-S[side][0])/(py-S[side][1])*(j-S[side][1]);
    I[1]=j;
    break;
  }
}

/* polyclip_area - Compute area of a given polygon (peform in double) */
float polyclip_area(float *px,float *py, int n) {
  int i,k;
  double area=0.0;
  for(i=0;i<n;i++) {
    k=(i==n-1?0:i+1);
    area+=(double)px[i]*(double)py[k]-(double)py[i]*(double)px[k];
  }
  if(area<0.) area=-area;
  return area/2.0;
}
