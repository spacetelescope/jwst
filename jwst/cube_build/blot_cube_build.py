""" Main module for blotting sky cube back to detector space
"""
import numpy as np
import logging
import time
from .. import datamodels
from ..assign_wcs import nirspec
from gwcs import wcstools
from . import instrument_defaults
from .blot_median import blot_wrapper  # c extension
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CubeBlot():

    def __init__(self, median_model, input_models):
        """Class Blot holds the main varibles for blotting sky cube to detector

        Information is pulled out of the median sky cube created by a previous
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
           The input models used to create the median sky cube.

        Returns
        -------
        CubeBlot class initialzied
        """

        # Pull out the needed information from the Median IFUCube
        self.median_skycube = median_model
        self.instrument = median_model.meta.instrument.name

        # basic information about the type of data
        self.grating = None
        self.filter = None
        self.subchannel = None
        self.channel = None
        self.par_median = None

        if self.instrument == 'MIRI':
            self.channel = median_model.meta.instrument.channel
            self.subchannel = median_model.meta.instrument.band.lower()
            self.par_median = self.channel

        elif self.instrument == 'NIRSPEC':
            self.grating = median_model.meta.instrument.grating
            self.filter = median_model.meta.instrument.filter
            self.par_median = self.grating
        # ________________________________________________________________
        # set up x,y,z of Median Cube
        # Median cube shoud have linear wavelength
        xcube, ycube, zcube = wcstools.grid_from_bounding_box(
            self.median_skycube.meta.wcs.bounding_box,
            step=(1, 1, 1))

        # using wcs of ifu cube determine ra,dec,lambda
        self.cube_ra, self.cube_dec, self.cube_wave = \
            self.median_skycube.meta.wcs(xcube + 1, ycube + 1, zcube + 1)

        # pull out flux from the median sky cube that matches with
        # cube_ra,dec,wave
        self.cube_flux = self.median_skycube.data

        # remove all the nan values
        valid1 = ~np.isnan(self.cube_ra)
        valid2 = ~np.isnan(self.cube_dec)
        good_data = np.where(valid1 & valid2)
        self.cube_ra = self.cube_ra[good_data]
        self.cube_dec = self.cube_dec[good_data]
        self.cube_wave = self.cube_wave[good_data]
        self.cube_flux = self.cube_flux[good_data]

        # initialize blotted images to be original input images valid for the Median Image
        # read channel (MIRI) or grating (NIRSpec) value that the Median Image covers
        # only use input models in this range
        self.input_models = []
        for model in input_models:
            if self.instrument == 'MIRI':
                par = model.meta.instrument.channel
            if self.instrument == 'NIRSPEC':
                par = model.meta.instrument.grating

            found = par.find(self.par_median)
            if found > -1:
                self.input_models.append(model)
    # **********************************************************************

    def blot_info(self):
        """ Prints the basic paramters of the blot image and median sky cube
        """
        log.info('Information on Blotting')
        log.info('Working with instrument %s', self.instrument)
        log.info('shape of sky cube %f %f %f',
                 self.median_skycube.data.shape[2],
                 self.median_skycube.data.shape[1],
                 self.median_skycube.data.shape[0])

        if self.instrument == 'MIRI':
            log.info('Channel %s', self.channel)
            log.info('Sub-channel %s', self.subchannel)

        elif self.instrument == 'NIRSPEC':
            log.info('Grating %s', self.grating)
            log.info('Filter %s', self.filter)

        log.info('Number of input models %i ', len(self.input_models))

    # ***********************************************************************

    def blot_images(self):
        if self.instrument == 'MIRI':
            blotmodels = self.blot_images_miri()
        elif self.instrument == 'NIRSPEC':
            blotmodels = self.blot_images_nirspec()
        return blotmodels
    # ************************************************************************

    def blot_images_miri(self):
        """ Core blotting routine for MIRI

        This is the main routine for blotting the MIRI median sky cube back to
        the detector space and creating a blotting image for each input model
        1. Loop over every data model to be blotted and find ra,dec,wavelength
           for every pixel in a valid slice on the detector.
        2. Loop over every input model and using the inverse (backwards) transform
           convert the median sky cube values ra, dec, lambda to the blotted
           x, y detector value (x_cube, y_cube).
        3. For each input model loop over the blotted x,y values and find
           the x, y detector values that fall within the roi region.
           The blotted flux  = the weighted flux, where the weight is based on
           distance between the center of the blotted pixel and the detector pixel.
        """
        blot_models = datamodels.ModelContainer()
        instrument_info = instrument_defaults.InstrumentInfo()

        for model in self.input_models:
            blot = model.copy()
            blot.err = None
            blot.dq = None
            xstart = 0
            # ___________________________________________________________________
            # For MIRI we only work on one channel at a time
            this_par1 = self.channel

            # get the detector values for this model
            xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
            ysize, xsize = model.data.shape
            ydet, xdet = np.mgrid[:ysize, :xsize]
            ydet = ydet.flatten()
            xdet = xdet.flatten()
            self.ycenter_grid, self.xcenter_grid = np.mgrid[0:ysize,
                                                            0:xsize]

            xsize2 = xend - xstart + 1
            xcenter = np.arange(xsize2) + xstart
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
            valid_channel = np.logical_and(x_cube >= xstart, x_cube <= xend)
            valid_flux = (flux_cube != 0)

            fuse = np.where(valid & valid_channel & valid_flux)
            x_cube = x_cube[fuse]
            y_cube = y_cube[fuse]
            flux_cube = flux_cube[fuse]

            log.info('Blotting back to %s', model.meta.filename)
            # ______________________________________________________________________________
            # blot_wrapper is a c extension that finds:
            # the overlapping median cube spaxels with the detector pixels
            # A median spaxel that falls withing the ROI of the center of the
            # detector pixel  is flagged as an overlapping pixel.

            xsize2 = xcenter.shape[0]
            roi_det = 1.0  # Just large enough that we don't get holes

            t0 = time.time()
            # set up c wrapper for blotting
            result = blot_wrapper(roi_det, xsize, ysize, xstart, xsize2,
                                  xcenter, ycenter,
                                  x_cube, y_cube, flux_cube)
            blot_flux, blot_weight = result
            igood = np.where(blot_weight > 0)
            blot_flux[igood] = blot_flux[igood] / blot_weight[igood]
            blot_flux = blot_flux.reshape((ysize, xsize))

            t1 = time.time()
            log.debug(f"Time to blot median image to input model =  {t1-t0:.1f}")

            blot.data = blot_flux
            blot_models.append(blot)
        return blot_models

    # ************************************************************************

    def blot_images_nirspec(self):
        """ Core blotting routine for NIRSPEC

        This is the main routine for blotting the NIRSPEC median sky cube back to
        the detector space and creating a blotting image for each input model.
        This routine was split from the MIRI routine because the blotting for NIRSpec
        needs to be done slice by slice and an error in the inverse mapping (sky to
        detector) mapped too many values back to the detector. This routine adds
        a check and first pulls out the min and max ra and dec values in the slice
        and only inverts the slice values back to the detector.
        For each data model loop over the 30 slices and find:
        a. the x,y bounding box of slice
        b. the ra, dec, lambda values for the x,y pixels in the slice
        c. from step b, determine the min and max ra and dec for slice values
        d. pull out the valid ra, dec and lambda values from the median sky cube that
           fall within the min and max ra,dec determined in step c
        e. invert the valid ra,dec, and lambda values for the slice determined in
           step d to the detector.
        f. blot the inverted x,y values to the detector plane. This steps finds the
           determines the overlap of the blotted x,y values with a regular grid setup
           in the detector plane which is the blotted image.

        """
        blot_models = datamodels.ModelContainer()

        for model in self.input_models:
            blot_ysize, blot_xsize = model.shape
            ntotal = blot_ysize * blot_xsize
            blot_flux = np.zeros(ntotal, dtype=np.float32)
            blot_weight = np.zeros(ntotal, dtype=np.float32)
            blot = model.copy()
            blot.err = None
            blot.dq = None

            # ___________________________________________________________________
            ycenter = np.arange(blot_ysize)
            xcenter = np.arange(blot_xsize)

            # for NIRSPEC wcs information accessed seperately for each slice
            nslices = 30
            log.info('Blotting  30 slices on NIRSPEC detector, this takes a little while')
            t0 = time.time()
            roi_det = 1.0  # Just large enough that we don't get holes

            x_total = []
            y_total = []
            flux_total = []
            for ii in range(nslices):
                # for each slice pull out the blotted values that actually fall on the slice region
                # use the bounding box of each slice to determine the slice limits
                slice_wcs = nirspec.nrs_wcs_set_input(model, ii)
                x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
                # using forward transform to limit the ra and dec values to
                # invert. The inverse transform for NIRSpec maps too many
                ra, dec, lam = slice_wcs(x, y)
                ramin = np.nanmin(ra)
                ramax = np.nanmax(ra)
                decmin = np.nanmin(dec)
                decmax = np.nanmax(dec)
                use1 = np.logical_and(self.cube_ra >= ramin, self.cube_ra <= ramax)
                use2 = np.logical_and(self.cube_dec >= decmin, self.cube_dec <= decmax)
                use = np.logical_and(use1, use2)

                ra_use = self.cube_ra[use]
                dec_use = self.cube_dec[use]
                wave_use = self.cube_wave[use]
                flux_use = self.cube_flux[use]

                # median cube ra,dec,wave -> x_slice, y_slice
                x_slice, y_slice = slice_wcs.invert(ra_use, dec_use, wave_use)
                x_slice = np.ndarray.flatten(x_slice)
                y_slice = np.ndarray.flatten(y_slice)
                flux_slice = np.ndarray.flatten(flux_use)

                # only use values what fall in bounding box of the slice
                xlimit, ylimit = slice_wcs.bounding_box
                xuse = np.logical_and(x_slice >= xlimit[0], x_slice <= xlimit[1])
                yuse = np.logical_and(y_slice >= ylimit[0], y_slice <= ylimit[1])
                fgood = (flux_slice != 0)  # edge cube pixel may = 0

                # fuse = np.logical_and(xuse, yuse)
                fuse = np.where(fgood & xuse & yuse)
                x_slice = x_slice[fuse]
                y_slice = y_slice[fuse]
                flux_slice = flux_slice[fuse]

                nslice = len(x_slice)
                for i in range(nslice):
                    x_total.append(x_slice[i])
                    y_total.append(y_slice[i])
                    flux_total.append(flux_slice[i])

            # end looping over the 30 slices
            t05 = time.time()
            # set up c wrapper for blotting
            xstart = 0
            xsize2 = blot_xsize
            x_total = np.array(x_total)
            y_total = np.array(y_total)
            flux_total = np.array(flux_total)

            result = blot_wrapper(roi_det, blot_xsize, blot_ysize, xstart, xsize2,
                                  xcenter, ycenter,
                                  x_total, y_total, flux_total)
            blot_flux_slice, blot_weight_slice = result
            blot_flux = blot_flux + blot_flux_slice
            blot_weight = blot_weight + blot_weight_slice

            result = None
            blot_weight_slice = None
            blot_flux_slice = None
            # done mapping median cube  to this input model
            t1 = time.time()
            log.debug(f"Time to blot median image (all 30 slices)  to input model =  {t1-t0:.1f}")
            log.debug(f"Time in blot_wrapper  =  {t1-t05:.1f}")
            igood = np.where(blot_weight > 0)
            blot_flux[igood] = blot_flux[igood] / blot_weight[igood]
            blot_flux = blot_flux.reshape((blot_ysize, blot_xsize))
            blot.data = blot_flux
            blot_models.append(blot)
        return blot_models
