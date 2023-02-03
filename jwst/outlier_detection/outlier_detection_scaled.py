"""Class definition for performing outlier detection with scaling."""

from copy import deepcopy
import numpy as np

from photutils.aperture import (aperture_photometry, CircularAperture, CircularAnnulus)
import astropy.units as u

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from ..resample import resample_utils
from ..tso_photometry.tso_photometry import tso_aperture_photometry
from .outlier_detection import OutlierDetection

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class OutlierDetectionScaled(OutlierDetection):
    """Class definition for applying scaled outlier detection.

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
        """Initialize class with input_models.

        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        reffiles : dict of `jwst.datamodels.DataModel`
            Dictionary of datamodels.  Keys are reffile_types.

        """
        super().__init__(self, input_models, reffiles=reffiles, **pars)

    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""
        self.build_suffix(**self.outlierpars)
        self._convert_inputs()

        pars = self.outlierpars
        save_intermediate_results = pars['save_intermediate_results']

        # Start by performing initial TSO Photometry on stack of DataModels
        # TODO:  need information about the actual source position in
        # TSO imaging mode (for all subarrays).
        # Meanwhile, this is a placeholder representing the geometric
        # center of the image.
        nints, ny, nx = self.inputs.data.shape
        xcenter = (ny - 1) / 2.
        ycenter = (ny - 1) / 2.

        # all radii are in pixel units
        if self.inputs.meta.instrument.pupil == 'WLP8':
            radius = 50
            radius_inner = 60
            radius_outer = 70
        else:
            radius = 3
            radius_inner = 4
            radius_outer = 5

        apertures = CircularAperture((xcenter, ycenter), r=radius)
        aperture_mask = apertures.to_mask(method='center')
        # This mask has 1 for mask region, 0 for outside of mask
        median_mask = aperture_mask.to_image((ny, nx))
        inv_median_mask = np.abs(median_mask - 1)
        # Perform photometry
        catalog = tso_aperture_photometry(self.inputs, xcenter, ycenter,
                                          radius, radius_inner,
                                          radius_outer)

        # Extract net photometry for the source
        # This will be the value used for scaling the median image within
        # the aperture region
        phot_values = catalog['net_aperture_sum']

        # Convert CubeModel into ModelContainer of 2-D DataModels
        for image in self.input_models:
            image.wht = resample_utils.build_driz_weight(
                image,
                weight_type='ivm',
                good_bits=pars['good_bits']
            )

        # Initialize intermediate products used in the outlier detection
        input_shape = self.input_models[0].data.shape
        median_model = datamodels.ImageModel(init=input_shape)
        median_model.meta = deepcopy(self.input_models[0].meta)
        base_filename = self.inputs.meta.filename
        median_model.meta.filename = self.make_output_path(
            basepath=base_filename, suffix='median'
        )

        # Perform median combination on set of drizzled mosaics
        median_model.data = self.create_median(self.input_models)
        aper2 = CircularAnnulus((xcenter, ycenter), r_in=radius_inner,
                                r_out=radius_outer)

        tbl1 = aperture_photometry(median_model.data, apertures,
                                   error=median_model.data * 0.0 + 1.0)
        tbl2 = aperture_photometry(median_model.data, aper2,
                                   error=median_model.data * 0.0 + 1.0)

        aperture_sum = u.Quantity(tbl1['aperture_sum'][0])
        annulus_sum = u.Quantity(tbl2['aperture_sum'][0])
        annulus_mean = annulus_sum / aper2.area
        aperture_bkg = annulus_mean * apertures.area
        median_phot_value = aperture_sum - aperture_bkg

        if save_intermediate_results:
            log.info("Writing out MEDIAN image to: {}".format(
                     median_model.meta.filename))
            median_model.save(median_model.meta.filename)

        # Scale the median image by the initial photometry (only in aperture)
        # to create equivalent of 'blot' images
        # Area outside of aperture in median will remain unchanged
        blot_models = ModelContainer()
        for i in range(nints):
            scale_factor = float(phot_values[i] / median_phot_value)
            scaled_image = datamodels.ImageModel(init=median_model.data.shape)
            scaled_image.meta = deepcopy(median_model.meta)
            scaled_data = (median_model.data * (scale_factor * median_mask) + (
                           median_model.data * inv_median_mask))
            scaled_image.data = scaled_data
            blot_models.append(scaled_image)

        if save_intermediate_results:
            log.info("Writing out Scaled Median images...")

            def make_output_path(ignored, idx=None):
                output_path = self.make_output_path(
                    basepath=base_filename, suffix='blot', idx=idx,
                    component_format='_{asn_id}_{idx}'
                )
                return output_path

            blot_models.save(make_output_path)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        self.detect_outliers(blot_models)

        # clean-up (just to be explicit about being finished
        # with these results)
        del median_model, blot_models
