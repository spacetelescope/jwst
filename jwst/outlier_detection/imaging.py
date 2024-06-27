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

from .utils import _convert_inputs, _detect_outliers, _remove_file, create_median, save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_models,
    save_intermediate_results,
    good_bits,
    maskpt,
    snr,
    scale,
    backg,
    resample_data,
    weight_type,
    pixfrac,
    kernel,
    fillval,
    allowed_memory,
    in_memory,
    asn_id,
    make_output_path,
):
    """Flag outlier pixels in DQ of input images."""
    input_models = _convert_inputs(input_models, good_bits, weight_type)

    if resample_data:
        # Start by creating resampled/mosaic images for
        # each group of exposures
        output_path = make_output_path(basepath=input_models[0].meta.filename,
                        suffix='')
        output_path = os.path.dirname(output_path)
        resamp = resample.ResampleData(
            input_models,
            output=output_path,
            single=True,
            blendheaders=False,
            wht_type=weight_type,
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            good_bits=good_bits,
            in_memory=in_memory,
            asn_id=asn_id,
            allowed_memory=allowed_memory,
        )
        drizzled_models = resamp.do_drizzle(input_models)

    else:
        # for non-dithered data, the resampled image is just the original image
        drizzled_models = input_models
        for i in range(len(input_models)):
            drizzled_models[i].wht = build_driz_weight(
                input_models[i],
                weight_type=weight_type,
                good_bits=good_bits)

    # Initialize intermediate products used in the outlier detection
    with datamodel_open(drizzled_models[0]) as dm0:
        median_model = datamodels.ImageModel(dm0.data.shape)
        median_model.update(dm0)
        median_model.meta.wcs = dm0.meta.wcs

    # Perform median combination on set of drizzled mosaics
    median_model.data = create_median(drizzled_models, maskpt)

    if save_intermediate_results:
        save_median(median_model, make_output_path, asn_id)
    else:
        # since we're not saving intermediate results if the drizzled models
        # were written to disk, remove them
        if not in_memory:
            for fn in drizzled_models._models:
                _remove_file(fn)

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    _detect_outliers(
        input_models,
        median_model,
        snr,
        scale,
        backg,
        resample_data,
    )

    # clean-up (just to be explicit about being finished with
    # these results)
    del median_model
