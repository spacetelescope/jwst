"""Class definition for performing outlier detection on spectra."""
from functools import partial

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

from ..resample import resample_spec, resample_utils
from .outlier_detection import OutlierDetection

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionSpec"]


class OutlierDetectionSpec(OutlierDetection):
    """Class definition for performing outlier detection on spectra.

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
      5. Perform statistical comparison between blotted image and original
         image to identify outliers.
      6. Updates input data model DQ arrays with mask of detected outliers.

    """

    default_suffix = 's2d'

    def __init__(self, input_models, reffiles=None, **pars):
        """Initialize class with input_models.

        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        reffiles : dict of `stdatamodels.jwst.datamodels.JwstDataModel`
            Dictionary of datamodels.  Keys are reffile_types.

        pars : dict, optional
            Optional user-specified parameters to modify how
            outlier_detection will operate.  Valid parameters include:
            - resample_suffix

        """
        OutlierDetection.__init__(self, input_models,
                                  reffiles=reffiles, **pars)

    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""
        self._convert_inputs()
        self.build_suffix(**self.outlierpars)

        pars = self.outlierpars
        save_intermediate_results = pars['save_intermediate_results']
        if pars['resample_data'] is True:
            # Start by creating resampled/mosaic images for
            #  each group of exposures
            resamp = resample_spec.ResampleSpecData(self.input_models, single=True,
                                                    blendheaders=False, **pars)
            drizzled_models = resamp.do_drizzle()
            if save_intermediate_results:
                for model in drizzled_models:
                    model.meta.filename = self.make_output_path(
                        basepath=model.meta.filename,
                        suffix=self.resample_suffix
                    )
                    log.info("Writing out resampled spectra...")
                    model.save(model.meta.filename)
        else:
            drizzled_models = self.input_models
            for i in range(len(self.input_models)):
                drizzled_models[i].wht = resample_utils.build_driz_weight(
                    self.input_models[i],
                    weight_type='ivm',
                    good_bits=pars['good_bits'])

        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(drizzled_models[0].data.shape)
        median_model.meta = drizzled_models[0].meta
        median_model.meta.filename = self.make_output_path(
            basepath=self.input_models[0].meta.filename,
            suffix='median'
        )

        # Perform median combination on set of drizzled mosaics
        # create_median should be called as a method from parent class
        median_model.data = self.create_median(drizzled_models)

        if save_intermediate_results:
            log.info("Writing out MEDIAN image to: {}".format(
                     median_model.meta.filename))
            median_model.save(median_model.meta.filename)

        if pars['resample_data'] is True:
            # Blot the median image back to recreate each input image specified
            # in the original input list/ASN/ModelContainer
            blot_models = self.blot_median(median_model)
            if save_intermediate_results:
                log.info("Writing out BLOT images...")
                blot_models.save(
                    partial(self.make_output_path, suffix='blot')
                )
        else:
            # Median image will serve as blot image
            blot_models = ModelContainer()
            for i in range(len(self.input_models)):
                blot_models.append(median_model)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        self.detect_outliers(blot_models)

        # clean-up (just to be explicit about being finished
        #  with these results)
        del median_model, blot_models
