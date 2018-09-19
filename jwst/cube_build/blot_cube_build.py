# Routines used for building cubes
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
        """
        Short Summary
        -------------
        Class Blot cube holds all the main varibles for blotting an IFU Cube back to
        the detector
        Information is pulled out of the Median Sky Cube created by a previous run
        of cube_build in single mode.

        Basic parameters of the instrument the  data is for is stored.
        The ra,dec, and wavelenth of the median sky cube is set up

        Parameters
        ---------
        median_model: median input sky cube created from a median stack of all the
           individual input_models mapped to the full IFU cube imprint on the sky
        input_models: data model
           a blotted image is created for each science image

        Returns
        -------
        CubeBlot class initialzied
        """
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
        self.cdelt3 = median_model.meta.wcsinfo.cdelt3 * 3600.0
#_______________________________________________________________________________
        xcube, ycube, zcube = wcstools.grid_from_bounding_box(
            self.median_skycube.meta.wcs.bounding_box,
            step=(1, 1, 1))

        self.cube_ra, self.cube_dec, self.cube_wave = self.median_skycube.meta.wcs(
            xcube,
            ycube,
            zcube)

        flux = self.median_skycube.data
        self.cube_flux = flux
#wavelength slices
        self.lam_centers = self.cube_wave[:, 0, 0]
# initialize blotted images to be original input images

        self.input_models = input_models

#********************************************************************************

    def blot_info(self):
        """
        Short Summary
        ------------
        Prints the basic paramters of the blot image and median sky cube
        """
        log.info('Information on Blotting')
        log.info('Working with instrument %s %s', self.instrument, self.detector)
        log.info('shape of sky cube %f %f %f', self.naxis1, self.naxis2, self.naxis3)

        log.info('Instrument %s ', self.instrument)
        if self.instrument == 'MIRI':
            log.info('Channel %s', self.channel)
            log.info('Sub-channel %s', self.subchannel)

        elif self.instrument == 'NIRSPEC':
            log.info('Grating %s', self.grating)
            log.info('Filter %s', self.filter)
        log.info('ROI size (spatial and wave) %f %f', self.rois, self.roiw)
        log.info('Number of input models %i ', len(self.input_models))

#********************************************************************************
    def blot_images(self):
        """
        Short Summary
        ------------
        Core blotting module
        Initialize blot_model = input_model
        1. Loop over every data model to be blotted and find ra,dec,wavelength
           for every slice pixel.
        2. Using WCS of input image convert ra,dec of input model to tangent
           plane values: xi,eta
        3. For the median sky cube convert the ra,dec of each x,y in this cube
           to xi,eta using the wcs of the imput image
        4. a. Loop over  every input_model valid IFU slice pixel and find the
            median image pixels that fall within the ROI of the center of pixel.
            The ROI parameters are the same as those used to construct the
            median sky cube and are stored in the meta data of the Median Cube.
            b. After all the overlapping median pixels have been found find the
            weighted flux using this overlapping pixels. The weighting is based
            on the distance between the detector pixel and the median flux pixel
            in tangent plane plane. Additional weighting parameters are read in
            from the median cube meta data.
            c. The blotted flux (blot.data) = the weighted flux determined in
            step b.

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
#________________________________________________________________________________
# From the x,y pixel for detector. For MIRI we only work on one channel at a time

            if self.instrument == 'MIRI':
                this_par1 = self.channel # only one channel is blotted at a time
                ch_name = '_ch' + this_par1
                blot.meta.filename = filename[:indx] + ch_name + '_blot.fits'

                # get the detector values for this model
                xstart, xend = instrument_info.GetMIRISliceEndPts(this_par1)
                ydet, xdet = np.mgrid[:1024, :1032]

                #mask out the side channel we aren not working on
                pixel_mask = np.full(model.shape, False, dtype=bool)
                pixel_mask[:, xstart:xend] = True
                ra_det, dec_det, lam_det = model.meta.wcs(xdet, ydet)

            elif self.instrument == 'NIRSPEC':
                blot.meta.filename = filename[:indx] + '_blot.fits'
                # initialize the ra,dec, and wavelength arrays
                # we will loop over slices and fill in values
                # the flag_det will be set when a slice pixel is filled in
                #   at the end we will use this flag to pull out valid data
                ra_det = np.zeros((2048, 2048))
                dec_det = np.zeros((2048, 2048))
                lam_det = np.zeros((2048, 2048))
                flag_det = np.zeros((2048, 2048))

                # for NIRSPEC each file has 30 slices
                # wcs information access seperately for each slice

                nslices = 30
                log.info('Looping over 30 slices on NIRSPEC detector, this takes a little while')
                for ii in range(nslices):
                    slice_wcs = nirspec.nrs_wcs_set_input(model, ii)
                    x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
                    ra, dec, lam = slice_wcs(x, y)

                    # the slices are curved on detector so a rectangular region
                    # returns NaNs
                    valid = ~np.isnan(lam)
                    ra = ra[valid]
                    dec = dec[valid]
                    lam = lam[valid]
                    x = x[valid]
                    y = y[valid]

                    xind = _toindex(x)
                    yind = _toindex(y)
                    xind = np.ndarray.flatten(xind)
                    yind = np.ndarray.flatten(yind)
                    ra = np.ndarray.flatten(ra)
                    dec = np.ndarray.flatten(dec)
                    lam = np.ndarray.flatten(lam)
                    ra_det[yind, xind] = ra
                    dec_det[yind, xind] = dec
                    lam_det[yind, xind] = lam
                    flag_det[yind, xind] = 1

# Done looping over slices
            log.info('Blotting back %s', model.meta.filename)

            if self.instrument == 'MIRI':
                valid3 = np.isfinite(lam_det)
                good_data = valid3 & pixel_mask
            elif self.instrument == 'NIRSPEC':
                good_data = np.where(flag_det == 1)

            y, x = good_data
            ra_blot = ra_det[good_data]
            dec_blot = dec_det[good_data]
            wave_blot = lam_det[good_data]
            
            crval1 = model.meta.wcsinfo.crval1
            crval2 = model.meta.wcsinfo.crval2

            # x,y detector pixels --> xi, eta
            xi_blot, eta_blot = coord.radec2std(crval1, crval2,
                                               ra_blot, dec_blot)

            # cube spaxel ra,dec values --> xi, eta
            xi_cube, eta_cube = coord.radec2std(crval1, crval2,
                                               self.cube_ra, self.cube_dec)
            nplane = self.naxis1 * self.naxis2
            self.xi_centers = np.reshape(xi_cube[0, :, :], nplane)
            self.eta_centers = np.reshape(eta_cube[0, :, :], nplane)

            num = ra_blot.size
#________________________________________________________________________________
# For every detector pixel find the overlapping median cube spaxels.
# A median spaxel that falls withing the ROI of the center of the detector pixel
# in the tangent plane is flagged as an overlapping pixel

            for ipt in range(0, num - 1):
                # xx,yy are the index value of the orginal detector frame -
                # blot image
                yy = y[ipt]
                xx = x[ipt]
                # find the cube values that fall withing ROI of detector xx,yy
                xdistance = (xi_blot[ipt] - self.xi_centers)
                ydistance = (eta_blot[ipt] - self.eta_centers)
                radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
                # indexr holds the index of the sky median spaxels that fall within
                # the spatial  ROI of xx,yy location
                indexr = np.where(radius <= self.rois)
                #indexz holds the index of the sky median spaxels that fall within
                # the spectral ROI of wave length assocation with xx,yy
                indexz = np.where(abs(self.lam_centers - wave_blot[ipt]) <= self.roiw)
                # Pull out the Cube spaxels falling with ROI regions

                wave_found = self.lam_centers[indexz]
                xi_found = self.xi_centers[indexr]
                eta_found = self.eta_centers[indexr]
#________________________________________________________________________________
                # form the arrays to be used calculated the weighting
                d1 = np.array(xi_found - xi_blot[ipt]) / self.cdelt1
                d2 = np.array(eta_found - eta_blot[ipt]) / self.cdelt2
                d3 = np.array(wave_found - wave_blot[ipt]) / self.cdelt3

                dxy = d1 * d1 + d2 * d2
                dxy_matrix = np.tile(dxy[np.newaxis].T, [1, d3.shape[0]])
                d3_matrix = np.tile(d3 * d3, [dxy_matrix.shape[0], 1])

                wdistance = dxy_matrix + d3_matrix
                weight_distance = np.power(np.sqrt(wdistance), self.weight_power)
                weight_distance[weight_distance < lower_limit] = lower_limit
                weight_distance = 1.0 / weight_distance

                # determine the spaxel xx_cube,yy_cube values of these spaxels in
                # the ROI so they can be used to pull out the flux of the median
                # sky cube.
                yy_cube = (indexr[0] / self.naxis1).astype(np.int)
                xx_cube = indexr[0] - yy_cube * self.naxis1
                scf = np.array([self.cube_flux[zz, yy_cube[ir], xx_cube[ir]]
                                for ir, rr in enumerate(indexr[0]) for zz in indexz[0]])
                scf = np.reshape(scf, weight_distance.shape)
                blot_flux[yy, xx] = np.sum(weight_distance * scf)
                blot_weight = np.sum(weight_distance)

                # check for blot_weight !=0
                if blot_weight == 0:
                    blot_flux[yy, xx] = 0
                else:
                    blot_flux[yy, xx] = blot_flux[yy, xx] / blot_weight
#________________________________________________________________________________
            blot.data = blot_flux
            blot_models.append(blot)
        t1 = time.time()
        log.info("Time Blot images = %.1f.s" % (t1 - t0,))
        return blot_models

#********************************************************************************
