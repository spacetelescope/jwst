
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
#from gwcs.utils import _domain_to_bounds
from ..associations import Association
from .. import datamodels
from ..assign_wcs import nirspec
from ..assign_wcs import pointing
from . import instrument_defaults

from . import file_table
from . import spaxel
from . import cube_overlap
from . import cube_cloud
from . import data_types

from gwcs import wcstools


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class CubeBlot(object):
# CubeBlot - holds all the important information for Blotting an IFU Cube back to detecto:

    def __init__(self, median_model,
                 input_models):
        
        self.number_files = len(input_models)
        self.median_skycube = median_model
        self.instrument = median_model.meta.instrument.name
        self.detector = median_model.meta.instrument.detector

        self.grating = None
        self.filter = None
        self.subchannel = None
        self.channel = None
        self.rois = median_model.meta.roi_spatial
        self.roiw = median_model.meta.roi_wave
        self.naxis1 = self.median_skycube.data.shape[2]
        self.naxis2 = self.median_skycube.data.shape[1]
        self.naxis3 =  self.median_skycube.data.shape[0]

        if(self.instrument == 'MIRI'):
            self.channel = median_model.meta.instrument.channel
            self.subchannel = median_model.meta.instrument.band
        elif(self.instrument == 'NIRSPEC'):
            self.grating = median_model.meta.instrument.grating
            self.filter = median_model.meta.instrument.filter
            
# initialize blotted images to be original input images

        self.blot_models = datamodels.ModelContainer() 

        for model in input_models:
            blot = model.copy()  
            blot.err = None
            blot.dq = None

            filename = model.meta.filename
            indx = filename.rfind('.fits')
            blot.meta.filename = filename[:indx] + '_blot.fits' #set output name
            self.blot_models.append(blot)

#        self.Cdelt1 = median_model.meta.wcsinfo.cdelt1
#        self.Cdelt2 = median_model.meta.wcsinfo.cdelt2
#        self.Cdelt3 = median_model.meta.wcsinfo.cdelt3
#        self.Crpix1 = median_model.meta.wcsinfo.crpix1
#        self.Crpix2 = median_model.meta.wcsinfo.crpix2
#        self.Crpix3 = median_model.meta.wcsinfo.crpix3
#        self.Crval1 = median_model.meta.wcsinfo.crval1
#        self.Crval2 = median_model.meta.wcsinfo.crval2
#        self.Crval3 = median_model.meta.wcsinfo.crval3
        

        xcube,ycube,zcube = wcstools.grid_from_bounding_box(self.median_skycube.meta.wcs.bounding_box,
                                                step=(1,1,1))
        
        cube_pos1,cube_pos2,cube_pos3 = self.median_skycube.meta.wcs(xcube,ycube,zcube)
        num = self.naxis1* self.naxis2 * self.naxis3

        self.cube_ra = np.reshape(cube_pos1,num)
        self.cube_dec = np.reshape(cube_pos2,num)
        self.cube_wave = np.reshape(cube_pos3,num)

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
        log.info('Number of input models %i ',len(self.blot_models))
        
#********************************************************************************

    def blot_images(self):

        instrument_info = instrument_defaults.InstrumentInfo()
        print('First 10 elements cube - ra  ',self.cube_ra[0:10])
        print('First 10 elements cube - dec%',self.cube_dec[0:10])
        print('First 10 elements cube - wave',self.cube_wave[0:10])        
        for model in self.blot_models:

            if(self.instrument =='MIRI'):
                this_par1 = self.channel # only one channel is blotted at a time
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
                y, x = np.mgrid[:1024, xstart:xend]
                y = np.reshape(y, y.size)
                x = np.reshape(x, x.size)

                ra,dec,lam = model.meta.wcs(x,y)


#            outsci = np.zeros(model.shape,dtype=np.float32)
#            xdet,ydet = model.meta.wcs.backward_transform(self.cube_ra,
#                                                          self.cube_dec,
#                                                          self.cube_wave)
#            print('transform world to detector',xdet.shape)
#            print('x',xdet[0:50])
#            print('y',ydet[0:50])



#********************************************************************************

    def create_spaxel(self):
        """
        Short Summary
        -------------
        # now you have the size of cube - create an instance for each spaxel
        # create an empty spaxel list - this will become a list of Spaxel classses

        Parameter
        ----------

        Returns
        -------
        list of classes contained in spaxel
        """
#________________________________________________________________________________


        total_num = self.naxis1*self.naxis2*self.naxis3

        if(self.interpolation == 'pointcloud'):
            for t in range(total_num):
                self.spaxel.append(spaxel.Spaxel())
        else:
            for t in range(total_num):
                self.spaxel.append(spaxel.SpaxelAB())

        return self.spaxel

