
# Routines used for building cubes
import sys
import time
import numpy as np
import math
import json
import logging

from astropy.io import fits
from astropy.modeling import models

from ..associations import Association
from .. import datamodels
from ..assign_wcs import nirspec
from ..assign_wcs import pointing
from gwcs import wcstools
from . import instrument_defaults
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class CubeBlot():
# CubeBlot - holds all the important information for Blotting an IFU Cube back to detecto:

    def __init__(self, median_model,
                 input_models):
        
        #Pull out the needed information from the Median IFUCube
        self.median_skycube = median_model
        self.instrument = median_model.meta.instrument.name
        self.detector = median_model.meta.instrument.detector

        #information on how the IFUCube was constructed 
        self.weight_power = median_model.meta.ifu.weight_power
        self.weighting = median_model.meta.ifu.weighting # if AREA ABORT
        self.rois = median_model.meta.ifu.roi_spatial
        self.roiw = median_model.meta.ifu.roi_wave

        #basic information about the type of data
        self.grating = None
        self.filter = None
        self.subchannel = None
        self.channel = None

        if(self.instrument == 'MIRI'):
            self.channel = median_model.meta.instrument.channel
            self.subchannel = median_model.meta.instrument.band
        elif(self.instrument == 'NIRSPEC'):
            self.grating = median_model.meta.instrument.grating
            self.filter = median_model.meta.instrument.filter

        # find the ra,dec,lambda of each element of the IFUCube
        self.naxis1 = self.median_skycube.data.shape[2]
        self.naxis2 = self.median_skycube.data.shape[1]
        self.naxis3 =  self.median_skycube.data.shape[0]
        self.cdelt1 = median_model.meta.wcsinfo.cdelt1*3600.0
        self.cdelt2 = median_model.meta.wcsinfo.cdelt2*3600.0
        self.cdelt3 = median_model.meta.wcsinfo.cdelt3*3600.0

        xcube,ycube,zcube = wcstools.grid_from_bounding_box(self.median_skycube.meta.wcs.bounding_box,
                                                step=(1,1,1))
        
        cube_pos1,cube_pos2,cube_pos3 = self.median_skycube.meta.wcs(xcube,ycube,zcube)
        num = self.naxis1* self.naxis2 * self.naxis3
        flux = self.median_skycube.data
        self.cube_ra = cube_pos1
        self.cube_dec = cube_pos2
        self.cube_wave = cube_pos3
        self.cube_flux = flux

        self.lam_centers = self.cube_wave[:,0,0]
        

# initialize blotted images to be original input images

        self.input_models = input_models


#********************************************************************************
# Print basic parameters of blot images and blot median

    def blot_info(self):
        log.info('Information on Blotting')
        log.info('Working with instrument %s %s',self.instrument,self.detector)        
        log.info('shape of sky cube %f %f %f',self.naxis1,self.naxis2,self.naxis3)

        log.info('Instrument %s ',self.instrument)
        if(self.instrument=='MIRI'):
            log.info('Channel %s',self.channel)
            log.info('Sub-channel %s',self.subchannel)

        elif(self.instrument=='NIRSPEC'):
            log.info('Grating %s',self.grating)
            log.info('Filter %s',self.filter)
        log.info('ROI size (spatial and wave) %f %f',self.rois,self.roiw)
        log.info('Number of input models %i ',len(self.input_models))
        
#********************************************************************************

    def blot_images(self):

        blot_models = datamodels.ModelContainer() 
        lower_limit = 0.01
        instrument_info = instrument_defaults.InstrumentInfo()

        for model in self.input_models:
            blot = model.copy()  
            blot.err = None
            blot.dq = None

            filename = model.meta.filename
            indx = filename.rfind('.fits')

            blot_flux = np.zeros(model.shape,dtype=np.float32)
            blot_weight = np.zeros(model.shape,dtype=np.float32)
            blot_iflux = np.zeros(model.shape,dtype=np.float32)            

# From the x,y pixel for detector. For MIRI we only work on one channel at a time
 
            if self.instrument =='MIRI':
                this_par1 = self.channel # only one channel is blotted at a time
                ch_name = '_ch' + this_par1
                blot.meta.filename = filename[:indx] +ch_name+ '_blot.fits' #set output name

                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
                ydet,xdet=np.mgrid[:1024,:1032]

                #mask out the side channel we aren not working on 
                pixel_mask = np.full(model.shape,False,dtype=bool)
                pixel_mask[:,xstart:xend] = True

            elif self.instrument =='NIRSPEC':
                blot.meta.filename = filename[:indx] + '_blot.fits' #set output name
                ydet,xdet=np.mgrid[:2048,:2048]

# determine the ra and dec values for each x,y pixel on the detector.

            log.info('Blotting back %s',model.meta.filename)
            ra_det,dec_det,wave_det = model.meta.wcs(xdet,ydet)

            valid1 = np.isfinite(ra_det)
            valid2 = np.isfinite(dec_det)
            valid3 = np.isfinite(wave_det)
            if self.instrument == 'MIRI':
                value = valid1 & valid2 & valid3 & pixel_mask
            elif elf.instrument == 'NIRSPEC':
                value = valid1 & valid2 & valid3

            good_data =  np.where( value == True)
            y,x = good_data

            ra_blot = ra_det[good_data]
            dec_blot = dec_det[good_data]
            wave_blot = wave_det[good_data]
                
            crval1 = model.meta.wcsinfo.crval1
            crval2 = model.meta.wcsinfo.crval2
            xi_blot,eta_blot = coord.radec2std(crval1, crval2,ra_blot,dec_blot) # for each x,y pixel

            xi_cube,eta_cube = coord.radec2std(crval1, crval2, # cube values
                                                   self.cube_ra,self.cube_dec) 

            nplane = self.naxis1 * self.naxis2
            self.xi_centers = np.reshape(xi_cube[0,:,:],nplane)
            self.eta_centers =np.reshape(eta_cube[0,:,:],nplane)
                
            num = ra_blot.size
#                print('Size of pixel detectors looping over',num)
            iprint = 0 
#________________________________________________________________________________  
               # loop over the valid pixels on the detector
            for ipt in range(0, num - 1):
                # xx,yy are the index value of the orginal detector frame -
                # blot image
                yy = y[ipt]
                xx = x[ipt]                    
                    # find the cube values that fall withing ROI of detector xx,yy
                xdistance = (xi_blot[ipt] - self.xi_centers)
                ydistance = (eta_blot[ipt] -self.eta_centers)

                radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
                
                indexr = np.where(radius  <=self.rois)
                indexz = np.where(abs(self.lam_centers - wave_blot[ipt]) <= self.roiw)

                wave_found = self.lam_centers[indexz]        # z Cube values falling in wavelength roi
                xi_found = self.xi_centers[indexr]   # x Cube values within radius 
                eta_found = self.eta_centers[indexr]  # y cube values with the radius
#________________________________________________________________________________  
#loop over the median cube pixels that fall within the ROI
# of the detector center

                for iz, zz in enumerate(indexz[0]):
                    istart = zz * nplane
                    for ir, rr in enumerate(indexr[0]):
                        yy_cube = int(rr/self.naxis1)
                        xx_cube = rr - yy_cube*self.naxis1
                            
                        d1 = (xi_found[ir]- xi_blot[ipt])/self.cdelt1
                        d2 = (eta_found[ir] - eta_blot[ipt])/self.cdelt2
                        d3 = (wave_found[iz] - wave_blot[ipt])/self.cdelt3
                            
                        weight_distance = math.sqrt(d1*d1 + d2*d2 + d3*d3)
                        weight_distance = math.pow(weight_distance,self.weight_power)

                        if weight_distance < lower_limit: weight_distance = lower_limit
                        weight_distance = 1.0 / weight_distance


                        blot_flux[yy,xx] = blot_flux[yy,xx]  + \
                            weight_distance*self.cube_flux[zz,yy_cube,xx_cube]
                        blot_weight[yy,xx] =blot_weight[yy,xx]  + weight_distance
                        blot_iflux[yy,xx] = blot_iflux[yy,xx] + 1
#________________________________________________________________________________  
                        
                    # determine finnal flux    
                    if(blot_iflux[yy,xx] > 0) :
                        blot_flux[yy,xx] = blot_flux[yy,xx]/blot_weight[yy,xx]
                    
                        
                blot.data = blot_flux
            blot_models.append(blot)
        return blot_models

#********************************************************************************
