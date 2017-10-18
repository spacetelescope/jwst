from __future__ import (division, print_function, unicode_literals,
    absolute_import)

import time
import numpy as np
from collections import OrderedDict

from photutils import aperture_photometry, CircularAperture, CircularAnnulus
import astropy.units as u
from astropy.table import QTable
from stsci.image import median
from stsci.tools import bitmask
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

from .. import datamodels
from ..resample import resample, gwcs_blot
from ..tso_photometry.tso_photometry import tso_aperture_photometry
from .outlier_detection import OutlierDetection, CRBIT

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DEFAULT_SUFFIX = 'i2d'

#class OutlierDetectionScaled(object):
class OutlierDetectionScaled(OutlierDetection):
    """
    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed by
    the various functions and then controls the operation of this process
    through all the steps used for the detection.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model and merges
         them with any user-provided values
      2. Performs photometry on single source from each input observations
      3. Creates a normalized median image from all grouped observations
      4. Perform statistical comparison between scaled,normalized median image
         and original image to identify outliers.
      5. Updates input data model DQ arrays with mask of detected outliers.

    """

    def __init__(self, input_models, reffiles=None, **pars):
        """
        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        reffiles : dict of `jwst.datamodels.DataModel`
            Dictionary of datamodels.  Keys are reffile_types.
        """
        OutlierDetection.__init__(self, input_models, reffiles=reffiles, **pars)


    def do_detection(self):
        """Flag outlier pixels in DQ of input images
        """
        pars = self.outlierpars
        save_intermediate_results = pars['save_intermediate_results']

        # Start by performing initial TSO Photometry on stack of DataModels
        # TODO:  need information about the actual source position in
        # TSO imaging mode (for all subarrays).
        # Meanwhile, this is a placeholder representing the geometric
        # center of the image.
        nints, ny, nx = self.input_models.data.shape
        xcenter = (ny - 1) / 2.
        ycenter = (ny - 1) / 2.

        # all radii are in pixel units
        if self.input_models.meta.instrument.pupil == 'WLP8':
            radius = 50
            radius_inner = 60
            radius_outer = 70
        else:
            radius = 3
            radius_inner = 4
            radius_outer = 5

        apertures = CircularAperture((xcenter, ycenter), r=radius)
        aperture_mask = apertures.to_mask(method='center')[0]
        # This mask has 1 for mask region, 0 for outside of mask
        median_mask = aperture_mask.to_image((ny, nx))
        inv_median_mask = np.abs(median_mask - 1)
        # Perform photometry
        catalog = tso_aperture_photometry(self.input_models, xcenter, ycenter,
                                          radius, radius_inner,
                                          radius_outer)

        # Extract net photometry for the source
        # This will be the value used for scaling the median image within
        # the aperture region
        phot_values = catalog['net_aperture_sum']

        # Convert CubeModel into ModelContainer of 2-D DataModels
        input_models = datamodels.ModelContainer()
        for i in range(self.input_models.data.shape[0]):
            image = datamodels.ImageModel(data=self.input_models.data[i],
                    err=self.input_models.err[i], dq=self.input_models.dq[i])
            image.meta = self.input_models.meta
            image.wht = resample.build_driz_weight(image, wht_type='exptime', good_bits=pars['good_bits'])
            input_models.append(image)

        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(init=input_models[0].data.shape)
        median_model.meta = input_models[0].meta
        base_filename = self.input_models.meta.filename
        median_model.meta.filename = '_'.join(base_filename.split('_')[:2] +
            ['median.fits'])

        # Perform median combination on set of drizzled mosaics
        median_model.data = self.create_median(input_models)
        aper2 = CircularAnnulus((xcenter, ycenter), r_in=radius_inner,
                            r_out=radius_outer)

        tbl1 = aperture_photometry(median_model.data, apertures,
                                   error=median_model.data * 0.0 + 1.0)
        tbl2 = aperture_photometry(median_model.data, aper2,
                                   error=median_model.data * 0.0 + 1.0)

        aperture_sum = u.Quantity(tbl1['aperture_sum'][0])
        annulus_sum = u.Quantity(tbl2['aperture_sum'][0])
        annulus_mean = annulus_sum / aper2.area()
        aperture_bkg = annulus_mean * apertures.area()
        median_phot_value = aperture_sum - aperture_bkg

        if save_intermediate_results:
            log.info("Writing out MEDIAN image to: {}".format(median_model.meta.filename))
            median_model.save(median_model.meta.filename)

        # Scale the median image by the initial photometry (only in aperture)
        # to create equivalent of 'blot' images
        # Area outside of aperture in median will remain unchanged
        blot_models = datamodels.ModelContainer()
        for i in range(self.input_models.data.shape[0]):
            scale_factor = float(phot_values[i] / median_phot_value)
            scaled_image = datamodels.ImageModel(init=median_model.data.shape)
            scaled_image.meta = median_model.meta
            scaled_data = median_model.data * (scale_factor * median_mask) + \
                    (median_model.data * inv_median_mask)
            scaled_image.data = scaled_data
            blot_models.append(scaled_image)

        if save_intermediate_results:
            log.info("Writing out Scaled Median images...")
            blot_models.save()

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        self.detect_outliers(blot_models)

        for i in range(self.input_models.data.shape[0]):
            self.input_models.dq[i] = input_models[i].dq

        # clean-up (just to be explicit about being finished with these results)
        del median_model, blot_models
