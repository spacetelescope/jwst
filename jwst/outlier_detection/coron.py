"""
Main class for performing outlier detection.

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


import logging
import os

import numpy as np

from stdatamodels.jwst import datamodels

from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight, build_mask

from .utils import create_cube_median, flag_cr_update_model
from ._fileio import remove_file, save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_model,
    save_intermediate_results,
    good_bits,
    maskpt,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    weight_type,
    asn_id,
    make_output_path,
):
    """Flag outlier pixels in DQ of input images."""
    if not isinstance(input_model, datamodels.JwstDataModel):
        input_model = datamodels.open(input_model)

    if not isinstance(input_model, datamodels.CubeModel):
        raise Exception(f"Input must be a CubeModel: {input_model}")

    # FIXME don't store this on the model
    # FIXME weight_type could now be used here. Similar to tso data coron
    # data was previously losing var_rnoise due to the conversion from a cube
    # to a ModelContainer (which makes the default ivm weight ignore var_rnoise).
    # Now that it's handled as a cube we could use the var_rnoise.
    input_model.wht = build_mask(input_model.dq, good_bits).astype(np.float32)

    # Perform median combination on set of drizzled mosaics
    median_data = create_cube_median(input_model, maskpt)

    if save_intermediate_results:
        # make a median model
        median_model = datamodels.ImageModel(median_data)
        median_model.update(input_model)
        median_model.meta.wcs = input_model.meta.wcs

        save_median(median_model, make_output_path, asn_id)
        del median_model

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    flag_cr_update_model(
        input_model,
        median_data,
        snr1,
        snr2,
        scale1,
        scale2,
        backg,
        False,
    )
    return input_model
