from __future__ import (division, print_function, unicode_literals,
    absolute_import)


import time
import numpy as np
from collections import OrderedDict

from stsci.image import median
from stsci.tools import bitmask
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

from .. import datamodels
from ..resample import resample, gwcs_blot
from .outlier_detection import OutlierDetection, CRBIT, create_median
from ..cube_build.cube_build_step import CubeBuildStep
from ..cube_build import blot_cube_build


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DEFAULT_SUFFIX = 's3d'
cube_build_config = 'cube_build.cfg'

class OutlierDetectionIFU(OutlierDetection):
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
      2. Resamples all input images into grouped observation mosaics.
      3. Creates a median image from all grouped observation mosaics.
      4. Blot median image to match each original input image.
      5. Perform statistical comparison between blotted image and original image
         to identify outliers.
      6. Updates input data model DQ arrays with mask of detected outliers.

    """

    def __init__(self, input_models, 
                    reffiles=None,  
                    **pars):
        """
        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        drizzled_models : list of objects
            ModelContainer containing drizzled grouped input images

        reffiles : dict of `jwst.datamodels.DataModel`
            Dictionary of datamodels.  Keys are reffile_types.
                        
            
        """
        OutlierDetection.__init__(self, input_models, reffiles=reffiles, **pars)

    def do_detection(self):
        """Flag outlier pixels in DQ of input images
        """
        pars = self.outlierpars
        save_intermediate_results = pars['save_intermediate_results']

        # Start by creating resampled/mosaic images for each group of exposures
        #for ch in [1,2]:

        single_IFUCube_result = CubeBuildStep.call(self.input_models,
                                                   config_file=cube_build_config,
                                                   channel=ch,
                                                   single='true')

        for model in single_IFUCube_result:
            model.meta.filename += self.resample_suffix
            if save_intermediate_results:
                log.info("Writing out resampled IFU cubes...")
                model.save(model.meta.filename)

        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(init=single_IFUCube_result[0].data.shape)
        median_model.meta = single_IFUCube_result[0].meta
        base_filename = self.input_models[0].meta.filename
        median_model.meta.filename = '_'.join(base_filename.split('_')[:2] +
            ['median.fits'])

        # Perform median combination on set of drizzled mosaics
        median_model.data = create_median(single_IFUCube_result, **pars)

        if save_intermediate_results:
            log.info("Writing out MEDIAN image to: {}".format(median_model.meta.filename))
            median_model.save(median_model.meta.filename)

        if pars['resample_data'] is True:
            # Blot the median image back to recreate each input image specified in
            # the original input list/ASN/ModelContainer
            blot_models = blot_median(median_model, self.input_models, **pars)
            if save_intermediate_results:
                for model in blot_models:
                    log.info("Writing out BLOT images...")
                    model.save(model.meta.filename)
        else:
            # Median image will serve as blot image
            blot_models = datamodels.ModelContainer()
            for i in range(len(self.input_models)):
                blot_models.append(median_model)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        detect_outliers(self.input_models, blot_models,
            self.reffiles, **self.outlierpars)

        # clean-up (just to be explicit about being finished with these results)
        del median_model, blot_models



