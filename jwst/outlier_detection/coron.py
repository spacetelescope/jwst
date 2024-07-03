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

from stdatamodels.jwst.datamodels.util import open as datamodel_open
from stdatamodels.jwst import datamodels

from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight

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
    # TODO only allow CubeModel

    # FIXME don't store this on the model
    # FIXME does input have a var_rnoise?
    input_model.wht = build_driz_weight(input_model, weight_type=weight_type, good_bits=good_bits)

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
