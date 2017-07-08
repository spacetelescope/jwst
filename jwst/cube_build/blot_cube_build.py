
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
from . import file_table
from . import instrument_defaults
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
        print('size of median skycube',self.median_skycube.data.shape)


        self.blot_models = datamodels.ModelContainer()
        for model in input_models:
            blot_median = model.copy()
            blot_median.err = None
            blot_median.dq = None

            filename = model.meta.filename
            indx = filename.rfind('.fits')
            blot_median.meta.filename = filename[:indx] + '_blot.fits'
            self.blot_models.append(blot_median)

        self.instrument = median_model.meta.instrument.name
        self.detector = median_model.meta.instrument.detector
        print('Working with instrument',self.instrument,self.detector)
        
        self.Cdelt1 = median_model.meta.wcsinfo.cdelt1
        self.Cdelt2 = median_model.meta.wcsinfo.cdelt2
        self.Cdelt3 = median_model.meta.wcsinfo.cdelt3
        self.Crpix1 = median_model.meta.wcsinfo.crpix1
        self.Crpix2 = median_model.meta.wcsinfo.crpix2
        self.Crpix3 = median_model.meta.wcsinfo.crpix3
        self.Crval1 = median_model.meta.wcsinfo.crval1
        self.Crval2 = median_model.meta.wcsinfo.crval2
        self.Crval3 = median_model.meta.wcsinfo.crval3

        self.naxis1 = self.median_skycube.data.shape[1]
        self.naxis2 = self.median_skycube.data.shape[0]
        self.naxis3 =  self.median_skycube.data.shape[2]

        print('shape of sky cube',self.naxis1,self.naxis2,self.naxis3)
        self.a_min = 0
        self.a_max = 0
        self.b_min = 0
        self.b_max = 0
        self.lambda_min = 0
        self.lambda_max = 0
        self.xcoord = None
        self.ycoord = None
        self.zcoord = None


        xi_center,eta_center = coord.radec2std(self.Crval1, self.Crval2,ra_ave,dec_ave)
        
        xi_min,eta_min = coord.radec2std(self.Crval1, self.Crval2,ra_min,dec_min)
        xi_max,eta_max = coord.radec2std(self.Crval1, self.Crval2,ra_max,dec_max)


#********************************************************************************
    def blot_images(self):

         
        
        for model in self.blot_models:
            
            outsci = np.zeros(model.shape,dtype=np.float32)
            worldtov23 = model.meta.wcs.get_transform("world","v2v3")
            v2ab_transform = model.meta.wcs.get_transform('v2v3','alpha_beta')


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

