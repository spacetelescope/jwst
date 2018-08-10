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
        self.weighting = median_model.meta.ifu.weighting 
        self.rois = median_model.meta.ifu.roi_spatial
        self.roiw = median_model.meta.ifu.roi_wave

        #basic information about the type of data
        self.grating = None
        self.filter = None
        self.subchannel = None
        self.channel = None

        if self.instrument == 'MIRI':
            self.channel = median_model.meta.instrument.channel
            self.subchannel = median_model.meta.instrument.band
        elif self.instrument == 'NIRSPEC':
            self.grating = median_model.meta.instrument.grating
            self.filter = median_model.meta.instrument.filter

        # find the ra,dec,lambda of each element of the IFUCube
        self.naxis1 = self.median_skycube.data.shape[2]
        self.naxis2 = self.median_skycube.data.shape[1]
        self.naxis3 =  self.median_skycube.data.shape[0]
        self.cdelt1 = median_model.meta.wcsinfo.cdelt1*3600.0
        self.cdelt2 = median_model.meta.wcsinfo.cdelt2*3600.0
        self.cdelt3 = median_model.meta.wcsinfo.cdelt3*3600.0
#________________________________________________________________________________
        xcube,ycube,zcube = wcstools.grid_from_bounding_box(self.median_skycube.meta.wcs.bounding_box,
                                                            step=(1,1,1))
        
        cube_pos1,cube_pos2,cube_pos3 = self.median_skycube.meta.wcs(xcube,
                                                                     ycube,
                                                                     zcube)
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
        if self.instrument=='MIRI':
            log.info('Channel %s',self.channel)
            log.info('Sub-channel %s',self.subchannel)

        elif self.instrument=='NIRSPEC':
            log.info('Grating %s',self.grating)
            log.info('Filter %s',self.filter)
        log.info('ROI size (spatial and wave) %f %f',self.rois,self.roiw)
        log.info('Number of input models %i ',len(self.input_models))
        
#********************************************************************************

    def blot_images(self):

        t0 = time.time()
        blot_models = datamodels.ModelContainer() 
        lower_limit = 0.01
        print_limit = 100000
        instrument_info = instrument_defaults.InstrumentInfo()

        for model in self.input_models:
            blot = model.copy()  
            blot.err = None
            blot.dq = None

            filename = model.meta.filename
            indx = filename.rfind('.fits')

            blot_flux = np.zeros(model.shape,dtype=np.float32)
#________________________________________________________________________________
# From the x,y pixel for detector. For MIRI we only work on one channel at a time
 
            if self.instrument =='MIRI':
                this_par1 = self.channel # only one channel is blotted at a time
                ch_name = '_ch' + this_par1
                blot.meta.filename = filename[:indx] +ch_name+ '_blot.fits' 

                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
                ydet,xdet=np.mgrid[:1024,:1032]

                #mask out the side channel we aren not working on 
                pixel_mask = np.full(model.shape,False,dtype=bool)
                pixel_mask[:,xstart:xend] = True
                ra_det,dec_det,lam_det = model.meta.wcs(xdet,ydet)

            elif self.instrument =='NIRSPEC':
                blot.meta.filename = filename[:indx] + '_blot.fits' 
                ra_det=np.zeros((2048,2048))
                dec_det=np.zeros((2048,2048))
                lam_det=np.zeros((2048,2048))
                flag_det = np.zeros((2048,2048))
                
                # for NIRSPEC each file has 30 slices - 
                # wcs information access seperately for each slice
                start_slice = 0
                end_slice = 29
                nslices = end_slice - start_slice + 1
                regions = list(range(start_slice, end_slice + 1))
#                log.info('Looping over 30 slices on NIRSPEC detector, this takes a little while')
                for ii in regions:
                    slice_wcs = nirspec.nrs_wcs_set_input(model, ii)
                    x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
                    ra, dec, lam = slice_wcs(x, y)

                    valid1 = np.isfinite(ra)
                    valid2 = np.isfinite(dec)
                    valid3 = np.isfinite(lam)
                    value = valid1 & valid2 & valid3
                    good_data =  np.where( value == True)
                    
                    ra = ra[good_data]
                    dec = dec[good_data]
                    lam = lam[good_data]
                    x = x[good_data]
                    y = y[good_data]
                    
                    xind = x.astype(np.int)
                    yind = y.astype(np.int)
                    xind = np.ndarray.flatten(xind)
                    yind = np.ndarray.flatten(yind)
                    ra = np.ndarray.flatten(ra)
                    dec = np.ndarray.flatten(dec)
                    lam = np.ndarray.flatten(lam)
                    ra_det[yind,xind] = ra 
                    dec_det[yind,xind] = dec
                    lam_det[yind,xind] = lam
                    flag_det[yind,xind] = ra*0.0 + 1
# Done looping over slices 
            log.info('Blotting back %s',model.meta.filename)

            if self.instrument == 'MIRI':
                valid1 = np.isfinite(ra_det)
                valid2 = np.isfinite(dec_det)
                valid3 = np.isfinite(lam_det)
                value = valid1 & valid2 & valid3 & pixel_mask
                good_data =  np.where( value == True)
            elif self.instrument == 'NIRSPEC':
                good_data = np.where(flag_det == 1)

            y,x = good_data
            ra_blot = ra_det[good_data]
            dec_blot = dec_det[good_data]
            wave_blot = lam_det[good_data]
                
            crval1 = model.meta.wcsinfo.crval1
            crval2 = model.meta.wcsinfo.crval2

            xi_blot,eta_blot = coord.radec2std(crval1, crval2,
                                               ra_blot,dec_blot) 

            xi_cube,eta_cube = coord.radec2std(crval1, crval2, 
                                               self.cube_ra,self.cube_dec) 
            
            nplane = self.naxis1 * self.naxis2
            self.xi_centers = np.reshape(xi_cube[0,:,:],nplane)
            self.eta_centers =np.reshape(eta_cube[0,:,:],nplane)

            num = ra_blot.size
            iprint = 0

#________________________________________________________________________________  
#loop over the median cube pixels that fall within the ROI
# of the detector center
            ii = 0 
            for ipt in range(0, num - 1):
#                if ii == print_limit: 
#                    log.info('On point %i out of %i',ipt,num)
#                ii = ii +1 
#                if ii > print_limit: ii = 0
                # xx,yy are the index value of the orginal detector frame -
                # blot image
                yy = y[ipt]
                xx = x[ipt]  
                
                # find the cube values that fall withing ROI of detector xx,yy
                xdistance = (xi_blot[ipt] - self.xi_centers)
                ydistance = (eta_blot[ipt] -self.eta_centers)
                radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
                indexr = np.where(radius  <=self.rois)
                indexz = np.where(abs(self.lam_centers- wave_blot[ipt]) <= self.roiw)
                zdistance = wave_blot[ipt] - self.lam_centers

                # Find the Cube spaxels falling with ROI regions
                wave_found = self.lam_centers[indexz] 
                xi_found = self.xi_centers[indexr]    
                eta_found = self.eta_centers[indexr]  
#________________________________________________________________________________  
                d1 = np.array(xi_found- xi_blot[ipt])/self.cdelt1
                d2 = np.array(eta_found - eta_blot[ipt])/self.cdelt2
                d3 = np.array(wave_found - wave_blot[ipt])/self.cdelt3

                dxy = d1*d1 + d2*d2
                dxy_matrix = np.tile(dxy[np.newaxis].T, [1, d3.shape[0]])
                d3_matrix = np.tile(d3*d3, [dxy_matrix.shape[0], 1])

                wdistance = dxy_matrix + d3_matrix
                weight_distance = np.power(np.sqrt(wdistance), self.weight_power)
                weight_distance[weight_distance < lower_limit] = lower_limit
                weight_distance = 1.0 / weight_distance

                yy_cube = (indexr[0]/self.naxis1).astype(np.int)
                xx_cube = indexr[0] - yy_cube*self.naxis1

                scf = np.array([self.cube_flux[zz, yy_cube[ir], xx_cube[ir]] 
                                for ir, rr in enumerate(indexr[0]) for zz in indexz[0]]) 
                scf = np.reshape(scf, weight_distance.shape)

                blot_flux[yy,xx] = np.sum(weight_distance * scf)
                blot_weight = np.sum(weight_distance)

                # check for blot_weight !=0  
                if blot_weight == 0: 
                    blot_flux[yy,xx] = 0
                else:
                    blot_flux[yy,xx] = blot_flux[yy,xx] / blot_weight
#________________________________________________________________________________  
                        
                blot.data = blot_flux
            blot_models.append(blot)
        t1 = time.time()
        log.info("Time Blot images = %.1f.s" % (t1 - t0,))
        return blot_models

#********************************************************************************
