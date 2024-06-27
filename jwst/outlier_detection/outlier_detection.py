"""Primary code for performing outlier detection on JWST observations."""

import logging
import warnings
import os

import numpy as np

from stdatamodels.jwst.datamodels.util import open as datamodel_open
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight

from .utils import _remove_file, create_median, detect_outliers

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetection"]


class OutlierDetection:
    """Main class for performing outlier detection.

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

    def __init__(self, input_models, **pars):
        """
        Initialize the class with input ModelContainers.

        Parameters
        ----------
        input_models : list of DataModels, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        pars : dict, optional
            Optional user-specified parameters to modify how outlier_detection
            will operate.  Valid parameters include:
            - resample_suffix

        """
        self.outlierpars = {}
        self.outlierpars.update(pars)

    def _convert_inputs(self, inputs):
        """Convert input into datamodel required for processing.

        This base class works on imaging data, and relies on use of the
        ModelContainer class as the format needed for processing. However,
        the input may not always be a ModelContainer object, so this method
        will convert the input to a ModelContainer object for processing.
        Additionally, sub-classes may redefine this to set up the input as
        whatever format the sub-class needs for processing.

        """
        bits = self.outlierpars['good_bits']
        if isinstance(inputs, ModelContainer):
            return inputs
        input_models = ModelContainer()
        num_inputs = inputs.data.shape[0]
        log.debug("Converting CubeModel to ModelContainer with {} images".
                  format(num_inputs))
        for i in range(inputs.data.shape[0]):
            image = datamodels.ImageModel(data=inputs.data[i],
                                          err=inputs.err[i],
                                          dq=inputs.dq[i])
            image.meta = inputs.meta
            image.wht = build_driz_weight(image,
                                          weight_type=self.outlierpars['weight_type'],
                                          good_bits=bits)
            input_models.append(image)
        return input_models

    def do_detection(self, inputs):
        """Flag outlier pixels in DQ of input images."""

        self.resample_suffix = "_outlier_i2d.fits"  # TODO factor this out

        input_models = self._convert_inputs(inputs)

        pars = self.outlierpars

        if pars['resample_data']:
            # Start by creating resampled/mosaic images for
            # each group of exposures
            output_path = pars["make_output_path"](basepath=input_models[0].meta.filename,
                            suffix='')
            output_path = os.path.dirname(output_path)
            resamp = resample.ResampleData(input_models, output=output_path, single=True,
                                           blendheaders=False, **pars)
            drizzled_models = resamp.do_drizzle(input_models)

        else:
            # for non-dithered data, the resampled image is just the original image
            drizzled_models = input_models
            for i in range(len(input_models)):
                drizzled_models[i].wht = build_driz_weight(
                    input_models[i],
                    weight_type=pars['weight_type'],
                    good_bits=pars['good_bits'])

        # Initialize intermediate products used in the outlier detection
        with datamodel_open(drizzled_models[0]) as dm0:
            median_model = datamodels.ImageModel(dm0.data.shape)
            median_model.update(dm0)
            median_model.meta.wcs = dm0.meta.wcs

        # Perform median combination on set of drizzled mosaics
        median_model.data = create_median(drizzled_models, self.outlierpars['maskpt'])

        if self.outlierpars['save_intermediate_results']:
            self.save_median(median_model)
        else:
            # since we're not saving intermediate results if the drizzled models
            # were written to disk, remove them
            if not self.outlierpars['in_memory']:
                for fn in drizzled_models._models:
                    _remove_file(fn)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        detect_outliers(
            input_models,
            median_model,
            self.outlierpars["snr"],
            self.outlierpars["scale"],
            self.outlierpars["backg"],
            self.outlierpars["resample_data"],
        )

        # clean-up (just to be explicit about being finished with
        # these results)
        del median_model

    def save_median(self, median_model):
        '''
        Save median if requested by user

        Parameters
        ----------
        median_model : ~jwst.datamodels.ImageModel
            The median ImageModel or CubeModel to save
        '''
        if self.outlierpars.get('asn_id', None) is None:
            suffix_to_remove = self.resample_suffix
        else:
            suffix_to_remove = f"_{self.outlierpars['asn_id']}{self.resample_suffix}"
        median_model_output_path = self.outlierpars["make_output_path"](
            basepath=median_model.meta.filename.replace(suffix_to_remove, '.fits'),
            suffix='median')
        median_model.save(median_model_output_path)
        log.info(f"Saved model in {median_model_output_path}")
