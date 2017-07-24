
# Routines used for building cubes
from __future__ import absolute_import, print_function
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

class CubeBlot(object):
# CubeBlot - holds all the important information for Blotting an IFU Cube back to detecto:

    def __init__(self, median_model,
                 input_models):
        
        #Pull out the needed information from the Median IFUCube
        self.median_skycube = median_model
        self.instrument = median_model.meta.instrument.name
        self.detector = median_model.meta.instrument.detector
        #information on how the IFUCube was constructed 
        self.weight_power = median_model.meta.weight_power
        self.weighting = median_model.meta.weighting # if AREA ABORT
        self.rois = median_model.meta.roi_spatial
        self.roiw = median_model.meta.roi_wave

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

        xcube,ycube,zcube = wcstools.grid_from_bounding_box(self.median_skycube.meta.wcs.bounding_box,
                                                step=(1,1,1))
        
        cube_pos1,cube_pos2,cube_pos3 = self.median_skycube.meta.wcs(xcube,ycube,zcube)
        num = self.naxis1* self.naxis2 * self.naxis3
        flux = self.median_skycube.data
        self.cube_ra = np.reshape(cube_pos1,num)
        self.cube_dec = np.reshape(cube_pos2,num)
        self.cube_wave = np.reshape(cube_pos3,num)
        self.cube_flux = np.reshape(flux,num)            
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
            blot.meta.filename = filename[:indx] + '_blot.fits' #set output name
            print('Blotting back model',model.meta.filename)

            if(self.instrument =='MIRI'):
                
                this_par1 = self.channel # only one channel is blotted at a time
                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
#                xdet,ydet = wcstools.grid_from_bounding_box(model.meta.wcs.bounding_box, step=(1,1))
                ydet, xdet = np.mgrid[:1024, :1032]                                    

                blot_flux = np.zeros(model.shape,dtype=np.float32)
                blot_weight = np.zeros(model.shape,dtype=np.float32)
                blot_iflux = np.zeros(model.shape,dtype=np.float32)

                pixel_mask = np.full(model.shape,False,dtype=bool)
                pixel_mask[:,xstart:xend] = True
                ra_det,dec_det,wave_det = model.meta.wcs(xdet,ydet)

#                print('xstart and xend',xstart,xend)
#                print('size of blot flux',blot_flux.size)                

                valid1 = np.isfinite(ra_det)
                valid2 = np.isfinite(dec_det)
                valid3 = np.isfinite(wave_det)
                value = valid1 & valid2 & valid3 & pixel_mask
                good_data =  np.where( value == True)
                y,x = good_data

#                print('size of value',value.size,value.shape)                
#                print(value[0,0:100])
#                print('x', x[0:100])
#                print(valid1[10,0:100])
#                print(pixel_mask[10,0:100])

#                print('length of index ',len(good_data[0]))                       
#                print('good data',good_data)
                
                ra_blot = ra_det[good_data]
                dec_blot = dec_det[good_data]
                wave_blot = wave_det[good_data]
                
                crval1 = model.meta.wcsinfo.crval1
                crval2 = model.meta.wcsinfo.crval2
                xi_blot,eta_blot = coord.radec2std(crval1, crval2,ra_blot,dec_blot) # for each x,y pixel

                xi_cube,eta_cube = coord.radec2std(crval1, crval2, # cube values
                                                   self.cube_ra,self.cube_dec) 

                print('size of xi_blot',xi_blot.shape)
                num = ra_blot.size
                print('Size of pixel detectors looping over/1024',num/1024)
                iprint = 0 
#________________________________________________________________________________  
                # loop over the valid pixels on the detector
                for ipt in range(0, num - 1):
                    
                    # xx,yy are the index value of the orginal detector frame -
                    # blot image
                    yy = y[ipt]
                    xx = x[ipt]                    
                    # find the cube values that fall withing ROI of detector xx,yy
                    xdistance = (xi_blot[ipt] - xi_cube)
                    ydistance = (eta_blot[ipt] - eta_cube)
                    radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
                
                    index =  np.where(np.logical_and( 
                            (radius  <=self.rois),  
                            (abs(wave_blot[ipt] - self.cube_wave) <= self.roiw)
                            ))
#                    print('for detector point',ipt,xi_blot[ipt],eta_blot[ipt],wave_blot[ipt])
#                    print(' original x,y index',x[ipt],y[ipt])
#                    print(' index, # and value ',len(index[0]),index)


                    wave_found = self.cube_wave[index]        # z Cube values falling in wavelength roi
                    xi_found = xi_cube[index]   # x Cube values within radius 
                    eta_found = eta_cube[index]  # y cube values with the radius
#                    print('wave found',wave_found)
#                    print('xi_found',xi_found)
#                    print('eta found',eta_found)

#________________________________________________________________________________  
                    #loop over the median cube pixels that fall within the ROI
                    # of the detector center
                    if(iprint == 0):
                        print ('on pixel ipt',ipt,len(index[0]))
                        
#                    print(xi_found.shape)     
                    d1 = xi_found- xi_blot[ipt]
                    d2 = eta_found - eta_blot[ipt]
                    d3 = wave_found - wave_blot[ipt]
                    d = d1*d1 + d2*d2 + d3*d3

                    for ii, ij in enumerate(index[0]):
                        #print(ipt,ii,ij,xx,yy)
                        #print(xi_found[ii],xi_blot[ipt],eta_found[ii],eta_blot[ipt])
#                        d1 = xi_found[ii] - xi_blot[ipt]
#                        d2 = eta_found[ii] - eta_blot[ipt]
#                        d3 = wave_found[ii] - wave_blot[ipt]
#                        weight_distance = math.sqrt(d1*d1 + d2*d2 + d3*d3)
#                        weight_distance = math.pow(weight_distance,self.weight_power)
#                        if weight_distance < lower_limit: weight_distance = lower_limit
#                        weight_distance = 1.0 / weight_distance
                        #print('index',index[0][j])


                        weight_distance = math.sqrt(d[ii])
                        weight_distance = math.pow(weight_distance,self.weight_power)
                        if weight_distance < lower_limit: weight_distance = lower_limit
                        weight_distance = 1.0 / weight_distance

                        blot_flux[yy,xx] = blot_flux[yy,xx]  + weight_distance*self.cube_flux[ij]
                        blot_weight[yy,xx] =blot_weight[yy,xx]  + weight_distance
                        blot_iflux[yy,xx] = blot_iflux[yy,xx] + 1
                        
                    # determine finnal flux    
                    if(blot_iflux[yy,xx] > 0) :
                        blot_flux[yy,xx] = blot_flux[yy,xx]/blot_weight[yy,xx]
                    
                    iprint = iprint+1
                    if(iprint == 2048):
                        iprint = 0 
                        
                blot.data = blot_flux
            blot_models.append(blot)
        return blot_models

#********************************************************************************
