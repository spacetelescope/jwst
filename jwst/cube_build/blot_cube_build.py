""" Main module for blotting sky cube back to detector space
"""
import time
import numpy as np
import logging

#from jwst.transforms.models import _toindex
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
        self.cube_flux = self.median_skycube.data

        # wavelength slices
        self.lam_centers = self.cube_wave[:, 0, 0]

        # initialize blotted images to be original input images
        self.input_models = input_models

    # ***********************************************************************

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

            filename = model.meta.filename
            indx = filename.rfind('.fits')

            blot_flux = np.zeros(model.shape, dtype=np.float32)
            # ___________________________________________________________________
            # For MIRI we only work on one channel at a time

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
                blot.meta.filename = filename[:indx] + '_blot.fits'
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
                # for NIRSPEC each file has 30 slices
                # wcs information accessed seperately for each slice

                nslices = 30
                log.info('Looping over 30 slices on NIRSPEC detector, this takes a little while')
                for ii in range(nslices):
                    slice_wcs = nirspec.nrs_wcs_set_input(model, ii)
                    x_slice, y_slice = slice_wcs.invert(self.cube_ra, self.cube_dec, self.cube_wave)
                    # there could be values that are nan because they do not fall on this slice.
                    x_slice = np.ndarray.flatten(x_slice)
                    y_slice = np.ndarray.flatten(y_slice)
                    flux_slice = np.ndarray.flatten(self.cube_flux)
                    # print('number of original cube values mapped',flux_slice.size)
                    valid = ~np.isnan(y_slice)
                    x_slice = x_slice[valid]
                    y_slice = y_slice[valid]
                    flux_slice = flux_slice[valid]
                    # print('number of valid cube values mapped',flux_slice.size)
                    if ii == 0:
                        x_cube = x_slice
                        y_cube = y_slice
                        flux_cube = flux_slice
                    else:
                        x_cube = np.concatenate((x_cube, x_slice), axis=0)
                        y_cube = np.concatenate((y_cube, y_slice), axis=0)
                        flux_cube = np.concatenate((flux_cube, flux_slice), axis=0)

            log.info('Blotting back %s', model.meta.filename)

            # ______________________________________________________________________________
            # For every detector pixel find the overlapping median cube spaxels.
            # A median spaxel that falls withing the ROI of the center of the detector
            # pixel in the tangent plane is flagged as an overlapping pixel.
            # the Regular grid is on the x,y detector

            blotf = self.blot_overlap_quick(xcenter, ycenter, x_cube, y_cube,
                                            flux_cube, blot_flux)
            blot.data = blotf
            blot_models.append(blot)
        t1 = time.time()
        log.info("Time Blot images = %.1f.s" % (t1 - t0,))
        return blot_models
    # ________________________________________________________________________________

    def blot_overlap(self, xdet, ydet, x_cube, y_cube, flux_cube, blot_flux):
        # blot_overlap finds to overlap between the blotted sky values (x_cube, y_cube)
        # and the detector pixels.
        # The looping is done over the x,y detector pixels (regular grid). The
        # distance to the x_cube, y_cube array is determiend  for each pixel to to find
        # the values in the roi of the detector pixel.
        # It was determiend the since x_cube and y_cube are large this method takes too
        # long. For now this routine was kept - maybe a techinque can be found to
        # speed it up.

        #t0 = time.time()
        # looping over detector pixels finding blotted values within roi (take a long time)
        num = xdet.size
        # print('Number of pixels on detector going to loop over',num)
        # print('Number of elements on sky mapped to detector going to search over',x_cube.size)

        roi_det = 0.5  # 1/2  a pixel
        for ipt in range(0, num - 1):

            # find the cube values that fall withing ROI of detector xdet, ydet
            xdistance = (x_cube - xdet[ipt])
            ydistance = (y_cube - ydet[ipt])
            radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
            # indexr holds the index of the sky median spaxels that fall
            # within the spatial  ROI of xx,yy location
            indexr = np.where(radius <= roi_det)
            num_roi = len(indexr[0])
            if num_roi > 0:
                x_found = x_cube[indexr]
                y_found = y_cube[indexr]
                # form the arrays to be used calculated the weighting
                d1 = np.array(x_found - xdet[ipt])
                d2 = np.array(y_found - ydet[ipt])

                dxy = d1 * d1 + d2 * d2
                dxy = np.sqrt(dxy)
                weight_distance = np.exp(-dxy)
                blot_flux[ydet[ipt], xdet[ipt]] = np.sum(weight_distance * flux_cube[indexr])
                blot_weight = np.sum(weight_distance)

                # check for blot_weight !=0
                if blot_weight == 0:
                    blot_flux[ydet[ipt], xdet[ipt]] = 0
                else:
                    blot_flux[ydet[ipt], xdet[ipt]] =  \
                        (blot_flux[ydet[ipt], xdet[ipt]] / blot_weight)

        #  after all looping
        # t1 = time.time()
        # print('Time to call blot_overlap',t1 - t0)
        return blot_flux
    # ________________________________________________________________________________

    def blot_overlap_quick(self, xcenter, ycenter, x_cube, y_cube,
                           flux_cube, blot_flux):

        # blot_overlap_quick finds to overlap between the blotted sky values (x_cube, y_cube)
        # and the detector pixels. This is faster than blot_overlap.
        # Looping is done over irregular arrays (x_cube, y_cube) and mapping to xcenter, ycenter
        # is quicker than looping over detector x,y and finding overlap with large array for
        # x_cube, y_cube
        # t0 = time.time()
        blot_weight = np.zeros_like(blot_flux)
        num = x_cube.size

        # print('blot_overlap_quick Number of elements on sky mapped to detector going to loop over',num)
        # TODO need to move this value (roi_det)  to reference file  or at least in spec ?
        roi_det = 0.5  # 1/2 a pixel
        for ipt in range(0, num - 1):
            # search xcenter and ycenter seperately. These arrays are smallsh.
            # xcenter size = naxis1 (naxis1/2 for MIRI)
            # ycenter size = naxis2
            xdistance = np.absolute(x_cube[ipt] - xcenter)
            ydistance = np.absolute(y_cube[ipt] - ycenter)
            index_x = np.where(xdistance <= roi_det)
            index_y = np.where(ydistance <= roi_det)

            if len(index_x[0]) > 0 and len(index_y[0]) > 0:
                for ix in index_x:
                    for iy in index_y:
                        # form the arrays to be used calculated the weighting
                        d1 = np.array(ix - x_cube[ipt])
                        d2 = np.array(iy - y_cube[ipt])
                        dxy = d1 * d1 + d2 * d2
                        dxy = np.sqrt(dxy)
                        weight_distance = np.exp(-dxy)
                        weighted_flux = weight_distance * flux_cube[ipt]
                        blot_flux[iy, ix] = blot_flux[iy, ix] + weighted_flux
                        blot_weight[iy, ix] = blot_weight[iy, ix] + weight_distance
        # done mapping blotted x,y (x_cube, y_cube) to detector
        good = blot_weight > 0
        blot_flux[good] = blot_flux[good] / blot_weight[good]

        # t1 = time.time()
        # print('Time to do quick blotting',t1- t0)
        return blot_flux
