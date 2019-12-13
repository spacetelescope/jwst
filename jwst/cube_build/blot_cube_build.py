""" Main module for blotting sky cube back to detector space
"""
import time
import numpy as np
import logging

from .. import datamodels
from ..assign_wcs import nirspec
from gwcs import wcstools
from . import instrument_defaults

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
        # self.weight_power = median_model.meta.ifu.weight_power
        # self.weighting = median_model.meta.ifu.weighting
        self.rois = median_model.meta.ifu.roi_spatial
        # self.roiw = median_model.meta.ifu.roi_wave

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
        self.cube_flux = self.median_skycube.data

        # wavelength slices
        self.lam_centers = self.cube_wave[:, 0, 0]

        # initialize blotted images to be original input images
        self.input_models = input_models

    # **********************************************************************

    def blot_info(self):
        """ Prints the basic paramters of the blot image and median sky cube
        """
        log.info('Information on Blotting')
        log.info('Working with instrument %s %s', self.instrument,
                 self.detector)
        log.info('shape of sky cube %f %f %f',
                 self.median_skycube.data.shape[2],
                 self.median_skycube.data.shape[1],
                 self.median_skycube.data.shape[0])

        log.info('Instrument %s ', self.instrument)
        if self.instrument == 'MIRI':
            log.info('Channel %s', self.channel)
            log.info('Sub-channel %s', self.subchannel)

        elif self.instrument == 'NIRSPEC':
            log.info('Grating %s', self.grating)
            log.info('Filter %s', self.filter)
        log.info('Spatial ROI size  %f', self.rois)
        log.info('Number of input models %i ', len(self.input_models))

    # ************************************************************************
    def blot_images(self):
        """ Core blotting routine

        This is the main routine for blotting the median sky image back to
        the detector space and creating a blotting image for each input model
        1. Loop over every data model to be blotted and find ra,dec,wavelength
           for every pixel in a valid slice on the detector.
        2. Using WCS of input  image convert the ra and dec of input model
           to then tangent plane values: xi,eta.
        3. Loop over every input model and using the inverse (backwards) transfor
           convert the median sky cube values ra, dec, lambda to the blotted
           x, y detector value (x_cibe, y_cube).

        4. For each input model loop over the blotted x,y values and find
            The x, y detector values that fall within the roi region. We loop
            over the blotted values instead of the detector x, y (the more
            obvious method) because of speed. We can break down the x, y detector
            pixels into a regular grid represented by an array of xcenters and
            ycenters. Because we are really intereseted finding all those blotted
            pixels that fall with the roi of a detector pixel we need to
            keep track of the blot flux and blot weight as we loop.
            c. The blotted flux  = the weighted flux determined in
            step b and stored in blot.data
        """
        t0 = time.time()
        blot_models = datamodels.ModelContainer()
        instrument_info = instrument_defaults.InstrumentInfo()

        for model in self.input_models:
            blot = model.copy()
            blot.err = None
            blot.dq = None

            blot_flux = np.zeros(model.shape, dtype=np.float32)
            # ___________________________________________________________________
            # For MIRI we only work on one channel at a time

            if self.instrument == 'MIRI':
                # only one channel is blotted at a time
                this_par1 = self.channel

                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)

                ysize, xsize = model.data.shape
                ydet, xdet = np.mgrid[:ysize, :xsize]
                ydet = ydet.flatten()
                xdet = xdet.flatten()

                self.ycenter_grid, self.xcenter_grid = np.mgrid[0:ysize,
                                                                0:xsize]

                xsize = xend - xstart + 1
                xcenter = np.arange(xsize) + xstart
                ycenter = np.arange(ysize)

                valid_channel = np.logical_and(xdet >= xstart, xdet <= xend)

                xdet = xdet[valid_channel]
                ydet = ydet[valid_channel]
                # cube spaxel ra,dec values --> x, y on detector
                x_cube, y_cube = model.meta.wcs.backward_transform(self.cube_ra,
                                                                   self.cube_dec,
                                                                   self.cube_wave)
                x_cube = np.ndarray.flatten(x_cube)
                y_cube = np.ndarray.flatten(y_cube)
                flux_cube = np.ndarray.flatten(self.cube_flux)
                valid = ~np.isnan(y_cube)
                x_cube = x_cube[valid]
                y_cube = y_cube[valid]
                flux_cube = flux_cube[valid]
                valid_channel = np.logical_and(x_cube >= xstart, x_cube <= xend)
                x_cube = x_cube[valid_channel]
                y_cube = y_cube[valid_channel]
                flux_cube = flux_cube[valid_channel]
            # ___________________________________________________________________
            # For NIRSPEC we need to work 1 slice at a time
            elif self.instrument == 'NIRSPEC':
                # initialize the x_cube,y_cube, and flux_cube arrays
                # we will loop over slices and fill in values
                # the flag_det will be set when a slice pixel is filled in
                #   at the end we will use this flag to pull out valid data
                ysize, xsize = model.data.shape
                x_cube = np.zeros((ysize, xsize))
                y_cube = np.zeros((ysize, xsize))
                flux_cube = np.zeros((ysize, xsize))

                ycenter = np.arange(ysize)
                xcenter = np.arange(xsize)
                self.ycenter_grid,self.xcenter_grid = np.mgrid[0:ysize,
                                                               0:xsize]

                # for NIRSPEC each file has 30 slices
                # wcs information accessed seperately for each slice

                nslices = 30
                log.info('Looping over 30 slices on NIRSPEC detector, this takes a little while')
                for ii in range(nslices):
                    # for each slice pull out the blotted values that actually fall on the slice region
                    # used the bounding box of each slice to determine the slice limits
                    slice_wcs = nirspec.nrs_wcs_set_input(model, ii)
                    x_slice, y_slice = slice_wcs.invert(self.cube_ra, self.cube_dec, self.cube_wave)
                    print('org min max',np.amin(x_slice),np.amax(x_slice),np.amin(y_slice),np.amax(y_slice))
                    x_slice = np.ndarray.flatten(x_slice)
                    y_slice = np.ndarray.flatten(y_slice)
                    flux_slice = np.ndarray.flatten(self.cube_flux)
                    print('number of original values',flux_slice.size)
                    xlimit, ylimit = slice_wcs.bounding_box
                    print('on slice ',ii, xlimit[0],xlimit[1],ylimit[0],ylimit[1])

                    xuse = np.logical_and(x_slice >= xlimit[0], x_slice <= xlimit[1])
                    yuse = np.logical_and(y_slice >= ylimit[0], y_slice <= ylimit[1])
                    fuse = np.logical_and(xuse, yuse)
                    x_slice = x_slice[fuse]
                    y_slice = y_slice[fuse]
                    print('min max of xslice',np.amin(x_slice),np.amax(x_slice))
                    print('min max of yslice',np.amin(y_slice),np.amax(y_slice))

                    flux_slice = flux_slice[fuse]
                    
                    print('number for slice',flux_slice.size)
                    if ii == 0:
                        x_cube = x_slice
                        y_cube = y_slice
                        flux_cube = flux_slice
                    else:
                        x_cube = np.concatenate((x_cube, x_slice), axis=0)
                        y_cube = np.concatenate((y_cube, y_slice), axis=0)
                        flux_cube = np.concatenate((flux_cube, flux_slice), axis=0)

            log.info('Blotting back to %s', model.meta.filename)

            # ______________________________________________________________________________
            # For every detector pixel find the overlapping median cube spaxels.
            # A median spaxel that falls withing the ROI of the center of the detector
            # pixel in the tangent plane is flagged as an overlapping pixel.
            # the Regular grid is on the x,y detector

            self.blot_overlap_quick(xcenter, ycenter, x_cube, y_cube,
                                            flux_cube, blot_flux)
            blot.data = blot_flux
            blot_models.append(blot)
        t1 = time.time()
        log.info("Time Blot images = %.1f.s" % (t1 - t0,))
        return blot_models
    # ________________________________________________________________________________

    def blot_overlap_quick(self, xcenter, ycenter, x_cube, y_cube,
                           flux_cube, blot_flux):

        # blot_overlap_quick finds to overlap between the blotted sky values (x_cube, y_cube)
        # and the detector pixels.
        # Looping is done over irregular arrays (x_cube, y_cube) and mapping to xcenter, ycenter
        # is quicker than looping over detector x,y and finding overlap with large array for
        # x_cube, y_cube
        # t0 = time.time()
        blot_weight = np.zeros_like(blot_flux)

        blot_ysize, blot_xsize = blot_flux.shape

        blot_flux = np.ndarray.flatten(blot_flux)
        blot_weight = np.ndarray.flatten(blot_weight)
        num = x_cube.size


        print('blot_overlap_quick Number of elements on sky mapped to detector going to loop over',num)
        # TODO need to move this value (roi_det)  to reference file  or at least in spec ?
        roi_det = 0.5  # 1/2 a pixel
        for ipt in range(num):
            # search xcenter and ycenter seperately. These arrays are smallsh.
            # xcenter size = naxis1 on detector (for MIRI only 1/2 array) 
            # ycenter size = naxis2 on detector
            xdistance = np.absolute(x_cube[ipt] - xcenter)
            ydistance = np.absolute(y_cube[ipt] - ycenter)

            index_x = np.where(xdistance <= roi_det)
            index_y = np.where(ydistance <= roi_det)

            if len(index_x[0]) > 0 and len(index_y[0]) > 0:

                d1 = np.array(x_cube[ipt] - xcenter[index_x])
                d2 = np.array(y_cube[ipt] - ycenter[index_y])
                dxy = (d1 * d1) + (d2 * d2)
                dxy = np.sqrt(dxy)
                weight_distance = np.exp(-dxy)
                weighted_flux = weight_distance * flux_cube[ipt]
                index2d = [ iy * blot_xsize + ix for iy in index_y[0] for ix in index_x[0]] 

                blot_flux[index2d] = blot_flux[index2d] + weighted_flux
                blot_weight[index2d] = blot_weight[index2d] + weight_distance

                if(ipt < 10) :
                    print('x y cube',x_cube[ipt], y_cube[ipt])
                    print('xcenters',xcenter[index_x])
                    print('ycenters',ycenter[index_y])
                    print('test',self.xcenter_grid[index2d], self.ycenter_grid[index2d])
                # form 
            #    for ix in index_x:
            #        for iy in index_y:
                        # form the arrays to be used calculated the weighting
            #            d1 = np.array(xcenter[ix] - x_cube[ipt])
            #            d2 = np.array(ycenter[iy] - y_cube[ipt])
            #            dxy = d1 * d1 + d2 * d2
            #            dxy = np.sqrt(dxy)
            #            weight_distance = np.exp(-dxy)
            #            weighted_flux = weight_distance * flux_cube[ipt]
            #            blot_flux[ycenter[iy], xcenter[ix]] = blot_flux[ycenter[iy], xcenter[ix]] + weighted_flux
            #            blot_weight[ycenter[iy], xcenter[ix]] = blot_weight[ycenter[iy], xcenter[ix]] + weight_distance


        # done mapping blotted x,y (x_cube, y_cube) to detector
        good = blot_weight > 0
        blot_flux[good] = blot_flux[good] / blot_weight[good]

        # t1 = time.time()
        # print('Time to do quick blotting',t1- t0)
        return blot_flux
