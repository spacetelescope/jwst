""" Main module for blotting sky cube back to detector space
"""
import time
import numpy as np
import logging

from jwst.transforms.models import _toindex
from .. import datamodels
from ..assign_wcs import nirspec
from gwcs import wcstools
from . import instrument_defaults
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CubeBlot():

    def __init__(self, median_model, input_models):
        """Class Blot holds the main varibles for blotting sky cube to detector


        Information is pulled out of the median sky Cube created by a previous
        run of cube_build in single mode and stored in the ClassBlot.These
        variables include the WCS of median sky cube, the weighting parameters
        used to create this median sky image and basic information of the input
        data (instrument, channel, band, grating or filter).

        Parameters
        ---------
        median_model: ifucube model
           The median input sky cube is created from a median stack of all the
           individual input_models mapped to the full IFU cube imprint on the
           sky.
        input_models: data model
           The input models used to create the median sky image.

        Returns
        -------
        CubeBlot class initialzied
        """

        # Pull out the needed information from the Median IFUCube
        self.median_skycube = median_model
        self.instrument = median_model.meta.instrument.name
        self.detector = median_model.meta.instrument.detector

        # information on how the IFUCube was constructed
        self.weight_power = median_model.meta.ifu.weight_power
        self.weighting = median_model.meta.ifu.weighting
        self.rois = median_model.meta.ifu.roi_spatial
        self.roiw = median_model.meta.ifu.roi_wave

        # basic information about the type of data
        self.grating = None
        self.filter = None
        self.subchannel = None
        self.channel = None

        if self.instrument == 'MIRI':
            self.channel = median_model.meta.instrument.channel
            self.subchannel = median_model.meta.instrument.band.lower()
        elif self.instrument == 'NIRSPEC':
            self.grating = median_model.meta.instrument.grating
            self.filter = median_model.meta.instrument.filter

        # find the ra,dec,lambda of each element of the IFUCube
        self.naxis1 = self.median_skycube.data.shape[2]
        self.naxis2 = self.median_skycube.data.shape[1]
        self.naxis3 = self.median_skycube.data.shape[0]
        self.cdelt1 = median_model.meta.wcsinfo.cdelt1 * 3600.0
        self.cdelt2 = median_model.meta.wcsinfo.cdelt2 * 3600.0
        self.cdelt3 = median_model.meta.wcsinfo.cdelt3
        self.crval1 = median_model.meta.wcsinfo.crval1
        self.crval2 = median_model.meta.wcsinfo.crval2
        # ________________________________________________________________
        # set up x,y,z of Median Cube 
        # Median cube shoud have linear wavelength
        xcube, ycube, zcube = wcstools.grid_from_bounding_box(
            self.median_skycube.meta.wcs.bounding_box,
            step=(1, 1, 1))

        # using wcs of ifu cube determine ra,dec,lambda
        self.cube_ra, self.cube_dec, self.cube_wave = \
            self.median_skycube.meta.wcs(xcube, ycube, zcube)

        # pull out flux from the median sky cube that matches with
        # cube_ra,dec,wave
        self.cube_flux =  self.median_skycube.data

        # wavelength slices
        self.lam_centers = self.cube_wave[:, 0, 0]

        # initialize blotted images to be original input images
        self.input_models = input_models

# *******************************************************************************

    def blot_info(self):
        """ Prints the basic paramters of the blot image and median sky cube
        """
        log.info('Information on Blotting')
        log.info('Working with instrument %s %s', self.instrument,
                 self.detector)
        log.info('shape of sky cube %f %f %f', self.naxis1, self.naxis2,
                 self.naxis3)

        log.info('Instrument %s ', self.instrument)
        if self.instrument == 'MIRI':
            log.info('Channel %s', self.channel)
            log.info('Sub-channel %s', self.subchannel)

        elif self.instrument == 'NIRSPEC':
            log.info('Grating %s', self.grating)
            log.info('Filter %s', self.filter)
        log.info('ROI size (spatial and wave) %f %f', self.rois, self.roiw)
        log.info('Number of input models %i ', len(self.input_models))

# *****************************************************************************
    def blot_images(self):
        """ Core blotting routine

        This is the main routine for blotting the median sky image back to
        the detector space and creating a blotting image for each input model
        1. Loop over every data model to be blotted and find ra,dec,wavelength
           for every pixel in a valid slice on the detector.
        2. Using WCS of input  image convert the ra and dec of input model
           to then tangent plane values: xi,eta.
        3. For the median sky cube convert the ra and dec of each x,y in this
           cube to xi,eta using the wcs of the input  image.
        4. a. Loop over  every input_model valid IFU slice pixel and find the
            median image pixels that fall within the ROI of the center of
            pixel.The ROI parameters are the same as those used to construct
            the median sky cube and are stored in the meta data of the median
            cube.
            b. After all the overlapping median pixels have been found find the
            weighted flux using these overlapping pixels. The weighting is
            based on the distance between the detector pixel and the median
            flux pixel in tangent plane plane. Additional weighting parameters
            are read in from the median cube meta data.
            c. The blotted flux  = the weighted flux determined in
            step b and stored in blot.data
        """
        t0 = time.time()
        blot_models = datamodels.ModelContainer()
        lower_limit = 0.01
        instrument_info = instrument_defaults.InstrumentInfo()

        for model in self.input_models:
            blot = model.copy()
            blot.err = None
            blot.dq = None

            filename = model.meta.filename
            indx = filename.rfind('.fits')

            blot_flux = np.zeros(model.shape, dtype=np.float32)
# ________________________________________________________________________________
# From the x,y pixel for detector. For MIRI we only work on one channel at a time

            if self.instrument == 'MIRI':
                # only one channel is blotted at a time
                this_par1 = self.channel
                ch_name = '_ch' + this_par1
                blot.meta.filename = filename[:indx] + ch_name + '_blot.fits'

                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
                ysize, xsize = model.data.shape
                ydet, xdet = np.mgrid[:ysize, :xsize]
                ydet = ydet.flatten()
                xdet = xdet.flatten()
                print('x y det shape',xdet.shape,ydet.shape)
                
                valid_channel = np.logical_and(xdet >= xstart, xdet<= xend)
                xdet = xdet[valid_channel]
                ydet = ydet[valid_channel]
                print('x y det shape',xdet.shape,ydet.shape)

                # cube spaxel ra,dec values --> x, y on detector
                x_cube, y_cube = model.meta.wcs.backward_transform(self.cube_ra,
                                                                   self.cube_dec,
                                                                   self.cube_wave)

                print('y_cube',y_cube.size)
                print(self.cube_ra.shape)
                print(self.cube_flux.shape)
                print(y_cube.shape)
#                for i in range(645):
#                    yslice = y_cube[i,:,:]
#                    xslice = x_cube[i,:,:]
#                    valid = np.isfinite(yslice)
#                    print(yslice[valid], xslice[valid])
                    
                x_cube = np.ndarray.flatten(x_cube)
                y_cube = np.ndarray.flatten(y_cube)
                flux_cube = np.ndarray.flatten(self.cube_flux)
                
                valid = ~np.isnan(y_cube)
                x_cube = x_cube[valid]
                y_cube = y_cube[valid]
                flux_cube = flux_cube[valid]
                valid_channel = np.logical_and(x_cube >= xstart, x_cube<= xend)
                
                x_cube = x_cube[valid_channel]
                y_cube = y_cube[valid_channel]
                flux_cube = flux_cube[valid_channel]

            elif self.instrument == 'NIRSPEC':
                blot.meta.filename = filename[:indx] + '_blot.fits'
                # initialize the x_cube,y_cube, and flux_cube arrays
                # we will loop over slices and fill in values
                # the flag_det will be set when a slice pixel is filled in
                #   at the end we will use this flag to pull out valid data
                ysize, xsize = model.data.shape
                x_cube = np.zeros((ysize, xsize))
                y_cube = np.zeros((ysize, xsize))
                flux_cube = np.zeros((ysize, xsize))

                # for NIRSPEC each file has 30 slices
                # wcs information access seperately for each slice

                nslices = 30
                log.info('Looping over 30 slices on NIRSPEC detector, this takes a little while')
                for ii in range(nslices):
                    slice_wcs = nirspec.nrs_wcs_set_input(model, ii)                    
                    x_slice,y_slice = slice_wcs.invert(self.cube_ra,self.cube_dec,self.cube_wave)
                    # many values will be nan because they do not fall on this slice. 
                    
                    x_slice = np.ndarray.flatten(x_slice)
                    y_slice = np.ndarray.flatten(y_slice)
                    flux_slice = np.ndarray.flatten(self.cube_flux)
                
                    valid = ~np.isnan(y_slice)
                    x_slice = x_slice[valid]
                    y_slice = y_slice[valid]
                    flux_slice = flux_slice[valid]
                    if ii == 0: 
                        x_cube = x_slice
                        y_cube = y_slice
                        flux_cube  = flux_slice
                    else:
                        x_cube  = np.concatenate(x_cube, x_slice)
                        y_cube  = np.concatenate(y_cube, y_slice)
                        flux_cube  = np.concatenate(flux_cube, flux_slice)


            log.info('Blotting back %s', model.meta.filename)
            num = x_cube.size
            print('Number of elements on detector going to loop over',num)
            # ______________________________________________________________________________
            # For every detector pixel find the overlapping median cube spaxels.
            # A median spaxel that falls withing the ROI of the center of the detector
            # pixel in the tangent plane is flagged as an overlapping pixel.
            # the Regular grid is on the x,y detector 
            #
            num = xdet.size 
            roi_det = 0.5 # 1/2  a pixel
            for ipt in range(0, num - 1):

                # find the cube values that fall withing ROI of detector xdet, ydet
                xdistance = (x_cube - xdet[ipt])
                ydistance = (y_cube - ydet[ipt])
                radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
                # indexr holds the index of the sky median spaxels that fall
                # within the spatial  ROI of xx,yy location
                indexr = np.where(radius <= roi_det)

                x_found = x_cube[indexr]
                y_found = y_cube[indexr]
                #print('size of x found',x_found.shape)
# ______________________________________________________________________________
                # form the arrays to be used calculated the weighting
                d1 = np.array(x_found - xdet[ipt]) 
                d2 = np.array(y_found - ydet[ipt]) 

                dxy = d1 * d1 + d2 * d2
                dxy = np.sqrt(dxy)
                weight_distance = np.exp(-dxy)
                #print('weight distance',weight_distance.shape,flux_cube[indexr].shape)

                blot_flux[ydet[ipt], xdet[ipt]] = np.sum(weight_distance * flux_cube[indexr])
                blot_weight = np.sum(weight_distance)

                # check for blot_weight !=0
                if blot_weight == 0:
                    blot_flux[ydet[ipt], xdet[ipt]] = 0
                else:
                    blot_flux[ydet[ipt], xdet[ipt]] =  \
                        (blot_flux[ydet[ipt], xdet[ipt]] / blot_weight)
# ________________________________________________________________________________
            blot.data = blot_flux
            blot_models.append(blot)
        t1 = time.time()
        log.info("Time Blot images = %.1f.s" % (t1 - t0,))
        return blot_models
