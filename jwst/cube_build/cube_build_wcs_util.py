# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json
import os

from astropy.io import fits
from ..associations import load_asn
from .. import datamodels
from ..assign_wcs import nirspec
from . import coord

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
# HELPER ROUTINES for CubeData class defined in cube_build.py
# these methods relate to wcs type procedures.  
# determine_scale

#********************************************************************************
def determine_scale(self):
#********************************************************************************
    """
    Short Summary
    -------------
    Determine the scale (sampling) in the 3 dimensions for the cube
    If the IFU cube covers more than 1 band - then use the rules to 
    define the Spatial and Wavelength sample size to use for the cube
    Current Rule: using the minimum

    Parameters
    ----------
    self.instrument_info holds the defaults scales for each channel/subchannel (MIRI)
    or Grating (NIRSPEC)

    Returns
    -------
    scale, holding the scale for the 3 dimensions of the cube/

    """
    scale = [0, 0, 0]
    if self.instrument == 'MIRI':
        number_bands = len(self.band_channel)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_bands):
            this_channel = self.band_channel[i]
            this_sub = self.band_subchannel[i]
            a_scale, b_scale, w_scale = self.instrument_info.GetScale(this_channel,this_sub)
            if a_scale < min_a:
                min_a = a_scale
            if b_scale < min_b:
                min_b = b_scale
            if w_scale < min_w:
                min_w = w_scale
        scale = [min_a, min_b, min_w]

    elif self.instrument == 'NIRSPEC':
        number_gratings = len(self.band_grating)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_gratings):
            this_gwa = self.band_grating[i]
            a_scale, b_scale, w_scale = self.instrument_info.GetScale(this_gwa)
            if a_scale < min_a:
                min_a = a_scale
            if b_scale < min_b:
                min_b = b_scale
            if w_scale < min_w:
                min_w = w_scale
        scale = [min_a, min_b, min_w]

#________________________________________________________________________________
# check and see if the user has set the scale or set by cfg. 

    a_scale = scale[0]
    if self.scale1 != 0.0:
        a_scale = self.scale1

    b_scale = scale[1]
    if self.scale2 != 0.0:
        b_scale = self.scale2

    w_scale = scale[2]
        # temp fix for large cubes - need to change to variable wavelength scale
    if self.scalew == 0 and self.num_bands > 6:   
        w_scale  = w_scale*2            
    if self.scalew == 0 and self.num_bands > 9:   
        w_scale  = w_scale*2            
    if self.scalew != 0.0:
        w_scale = self.scalew

    scale = [a_scale, b_scale, w_scale]
    return scale
#_______________________________________________________________________


#********************************************************************************
def find_footprint_MIRI(self, input, this_channel, instrument_info):
#********************************************************************************

    """
    Short Summary
    -------------
    For each channel find:
    a. the min and max spatial coordinates (alpha,beta) or (V2-v3) depending on coordinate system.
      axis a = naxis 1, axis b = naxis2
    b. min and max wavelength is also determined. , beta and lambda for those slices


    Parameters
    ----------
    input: input model (or file)
    this_channel: channel working with


    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.
    spaxial coordinates are in units of arc seconds.
    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers

    xstart, xend = instrument_info.GetMIRISliceEndPts(this_channel)
    y, x = np.mgrid[:1024, xstart:xend]
    
    coord1 = np.zeros(y.shape)
    coord2 = np.zeros(y.shape)
    lam = np.zeros(y.shape)

    if self.coord_system == 'alpha-beta':
        detector2alpha_beta = input.meta.wcs.get_transform('detector', 'alpha_beta')

        shift = models.Shift(1) & models.Shift(1)
        for key in detector2alpha_beta.selector:
            detector2alpha_beta.selector[key] = shift | detector2alpha_beta.selector[key]

        input.meta.wcs.set_transform('detector','alpha_beta',detector2alpha_beta)

        coord1, coord2, lam = detector2alpha_beta(x, y)

#        xtest = 28.310396-1 # test pixel to compare with Distortion doc 
#        ytest = 512.0-1     # test pixel to compare with Distortion doc
#        coord1_test,coord2_test,lam_test = detector2alpha_beta(xtest,ytest)
#        print('test values',xtest+1,ytest+1,coord1_test,coord2_test,lam_test)


    elif self.coord_system == 'ra-dec':
        detector2v23 = input.meta.wcs.get_transform('detector', 'v2v3')
        v23toworld = input.meta.wcs.get_transform("v2v3","world")

        v2, v3, lam = detector2v23(x, y)
        coord1,coord2,lam = v23toworld(v2,v3,lam)

    else:
        # error the coordinate system is not defined
        raise NoCoordSystem(" The output cube coordinate system is not definded")

    a_min = np.nanmin(coord1)
    a_max = np.nanmax(coord1)

    b_min = np.nanmin(coord2)
    b_max = np.nanmax(coord2)

    lambda_min = np.nanmin(lam)
    lambda_max = np.nanmax(lam)

    return a_min, a_max, b_min, b_max, lambda_min, lambda_max

#********************************************************************************
def find_footprint_NIRSPEC(self, input,flag_data):
#********************************************************************************

    """
    Short Summary
    -------------
    For each slice find:
    a. the min and max spatial coordinates (alpha,beta) or (V2-v3) depending on coordinate system.
      axis a = naxis 1, axis b = naxis2
    b. min and max wavelength is also determined. , beta and lambda for those slices


    Parameters
    ----------
    input: input model (or file)

    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.

    """
    # loop over all the region (Slices) in the Channel
    # based on regions mask (indexed by slice number) find all the detector
    # x,y values for slice. Then convert the x,y values to  v2,v3,lambda
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers

    start_slice = 0
    end_slice = 29

    nslices = end_slice - start_slice + 1

    a_slice = np.zeros(nslices * 2)
    b_slice = np.zeros(nslices * 2)
    lambda_slice = np.zeros(nslices * 2)

    regions = list(range(start_slice, end_slice + 1))
    k = 0

    self.log.info('Looping over slices to determine cube size .. this takes a while')
    # for NIRSPEC there are 30 regions
    for i in regions:

        slice_wcs = nirspec.nrs_wcs_set_input(input,  i)
        yrange_slice = slice_wcs.bounding_box[1][0],slice_wcs.bounding_box[1][0]
        xrange_slice = slice_wcs.bounding_box[0][0],slice_wcs.bounding_box[0][0]

        if(xrange_slice[0] >= 0 and xrange_slice[1] > 0):

            x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box,step=(1,1), center=True)
            ra,dec,lam = slice_wcs(x,y)

            a_slice[k] = np.nanmin(ra)
            a_slice[k + 1] = np.nanmax(ra)

            b_slice[k] = np.nanmin(dec)
            b_slice[k + 1] = np.nanmax(dec)

            lambda_slice[k] = np.nanmin(lam)
            lambda_slice[k + 1] = np.nanmax(lam)

        k = k + 2

    a_min = min(a_slice)
    a_max = max(a_slice)

    b_min = min(b_slice)
    b_max = max(b_slice)

    lambda_min = min(lambda_slice)
    lambda_max = max(lambda_slice)

    if(a_min == 0.0 and a_max == 0.0 and b_min ==0.0 and b_max == 0.0):
        self.log.info('This NIRSPEC exposure has no IFU data on it - skipping file')
        flag_data = -1

    return a_min, a_max, b_min, b_max, lambda_min, lambda_max

#_______________________________________________________________________
# Footprint values are RA,DEC values on the sky
# Values are given in degrees 

def set_geometry(self, footprint):

        deg2rad = math.pi/180.0
        ra_min, ra_max, dec_min, dec_max,lambda_min, lambda_max = footprint # in degrees
        dec_ave = (dec_min + dec_max)/2.0

        # actually this is hard due to converenge of hour angle
        # improve determining ra_ave in the future - do not just average (BAD) 
        ra_ave = ((ra_min + ra_max)/2.0 )#* math.cos(dec_ave*deg2rad) 
#        range_ra = (ra_max - ra_min) * 3600.0 * math.cos(dec_ave*deg2rad)
#        range_dec = (dec_max - dec_min) * 3600.0

        self.Crval1 = ra_ave 
        self.Crval2 = dec_ave
        xi_center,eta_center = coord.radec2std(self.Crval1, self.Crval2,ra_ave,dec_ave)
        
        xi_min,eta_min = coord.radec2std(self.Crval1, self.Crval2,ra_min,dec_min)
        xi_max,eta_max = coord.radec2std(self.Crval1, self.Crval2,ra_max,dec_max)
#________________________________________________________________________________
        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        # to find location of center abs of min values is how many pixels 

        n1a = int(math.ceil(math.fabs(xi_min) / self.Cdelt1)) 
        n2a = int(math.ceil(math.fabs(eta_min) / self.Cdelt2)) 

        n1b = int(math.ceil(math.fabs(xi_max) / self.Cdelt1)) 
        n2b = int(math.ceil(math.fabs(eta_max) / self.Cdelt2)) 

        xi_min = 0.0 - (n1a * self.Cdelt1) - self.Cdelt1/2.0
        xi_max = (n1b * self.Cdelt1) + self.Cdelt1/2.0

        eta_min = 0.0 - (n2a * self.Cdelt2) - self.Cdelt2/2.0
        eta_max = (n2b * self.Cdelt2) + self.Cdelt2/2.0
        
        self.Crpix1 = n1a
        self.Crpix2 = n2a

        self.naxis1 = n1a + n1b
        self.naxis2 = n2a + n2b

        self.a_min  = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max

# center of spaxels 
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.Cdelt1 / 2.0
        for i in range(self.naxis1):
            self.xcoord[i] = xstart
            xstart = xstart + self.Cdelt1


        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.Cdelt2 / 2.0

        for i in range(self.naxis2):
            self.ycoord[i] = ystart
            ystart = ystart + self.Cdelt2
#            print('ycoord',self.ycoord[i],i)

#        print('ycoord xcoord shape',self.ycoord.shape,self.xcoord.shape)
        
#_______________________________________________________________________
        
#        ystart = self.ycoord[0]
#        yend = self.ycoord[0] + self.Cdelt2*(self.naxis2)
        
#        xstart = self.xcoord[0]
#        xend = self.xcoord[0] + self.Cdelt1*(self.naxis1)

#        yy,xx = np.mgrid[ystart:yend:self.Cdelt2,
#                         xstart:xend:self.Cdelt1]

        ygrid = np.zeros(self.naxis2*self.naxis1)
        xgrid = np.zeros(self.naxis2*self.naxis1)

        k = 0 
        ystart = self.ycoord[0]
        for i in range(self.naxis2):
            xstart = self.xcoord[0]
            for j in range(self.naxis1): 
                xgrid[k] = xstart
                ygrid[k] = ystart
                xstart = xstart + self.Cdelt1
                k = k + 1
            ystart = ystart + self.Cdelt2
        

#        print('y start end',ystart,yend)
#        print('x start end',xstart,xend)

#        print('yy shape',yy.shape,self.ycoord.shape)
#        print('xx shape',xx.shape,self.xcoord.shape)

#        self.Ycenters = np.ravel(yy)
#        self.Xcenters = np.ravel(xx)

        self.Xcenters = xgrid
        self.Ycenters = ygrid
#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        range_lambda = self.lambda_max - self.lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
        lambda_center = (self.lambda_max + self.lambda_min) / 2.0
        self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
        self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

        self.zcoord = np.zeros(self.naxis3)
        self.Crval3 = self.lambda_min
        self.Crpix3 = 1.0
        zstart = self.lambda_min + self.Cdelt3 / 2.0

        for i in range(self.naxis3):
            self.zcoord[i] = zstart
            zstart = zstart + self.Cdelt3


#_______________________________________________________________________
# cube in alpha-beta space (single exposure cube - small FOV assume rectangular coord system
def set_geometryAB(self, footprint):
    self.a_min, self.a_max, self.b_min, self.b_max, self.lambda_min, self.lambda_max = footprint

        #set up the a (x) coordinates of the cube
    range_a = self.a_max - self.a_min
    self.naxis1 = int(math.ceil(range_a / self.Cdelt1))

        # adjust min and max based on integer value of naxis1
    a_center = (self.a_max + self.a_min) / 2.0
    self.a_min = a_center - (self.naxis1 / 2.0) * self.Cdelt1
    self.a_max = a_center + (self.naxis1 / 2.0) * self.Cdelt1

    self.xcoord = np.zeros(self.naxis1)
    self.Crval1 = self.a_min
    self.Crpix1 = 0.5
    xstart = self.a_min + self.Cdelt1 / 2.0
    for i in range(self.naxis1):
        self.xcoord[i] = xstart
        xstart = xstart + self.Cdelt1

#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

    range_lambda = self.lambda_max - self.lambda_min
    self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
    lambda_center = (self.lambda_max + self.lambda_min) / 2.0

    self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
    self.lambda_max = lambda_center + (self.naxis3 / 2.0) * self.Cdelt3

    self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

    self.zcoord = np.zeros(self.naxis3)
    self.Crval3 = self.lambda_min
    self.Crpix3 = 1.0
    zstart = self.lambda_min + self.Cdelt3 / 2.0

    for i in range(self.naxis3):
        self.zcoord[i] = zstart
        zstart = zstart + self.Cdelt3
#_______________________________________________________________________

    range_b = self.b_max - self.b_min

    self.naxis2 = int(math.ceil(range_b / self.Cdelt2))
    b_center = (self.b_max + self.b_min) / 2.0

        # adjust min and max based on integer value of naxis2
    self.b_max = b_center + (self.naxis2 / 2.0) * self.Cdelt2
    self.b_min = b_center - (self.naxis2 / 2.0) * self.Cdelt2


    self.ycoord = np.zeros(self.naxis2)
    self.Crval2 = self.b_min
    self.Crpix2 = 0.5
    ystart = self.b_min + self.Cdelt2 / 2.0
    for i in range(self.naxis2):
        self.ycoord[i] = ystart
        ystart = ystart + self.Cdelt2



#_______________________________________________________________________
def print_cube_geometry(self):
        log.info('Cube Geometry:')
        blank = '  '
        if (self.coord_system == 'alpha-beta'):
            log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (alpha,beta arc sec)')
        else:
            log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (xi,eta arc sec)')
        log.info('Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis1, self.Crpix1, self.Crval1, self.Cdelt1, self.a_min, self.a_max)
        log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis2, self.Crpix2, self.Crval2, self.Cdelt2, self.b_min, self.b_max)
        log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis3, self.Crpix3, self.Crval3, self.Cdelt3, self.lambda_min, self.lambda_max)

        if(self.instrument == 'MIRI'):
            # length of channel and subchannel are the same 
            number_bands = len(self.band_channel)

            for i in range(number_bands):
                this_channel = self.band_channel[i]
                this_subchannel = self.band_subchannel[i]
                log.info('Cube covers channel, subchannel: %s %s ', this_channel,this_subchannel)
        elif(self.instrument == 'NIRSPEC'):
            # number of filters and gratings are the same
            number_bands = len(self.band_filter)

            for i in range(number_bands):
                this_fwa = self.band_filter[i]
                this_gwa = self.band_grating[i]
                log.info('Cube covers grating, filter: %s %s ', this_gwa,this_fwa)
