"""
Submodule for performing outlier detection on imaging data.
"""

import copy
import logging
import os

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary
from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight
from jwst.stpipe.utilities import record_step_status

from .utils import create_median_library, flag_crs_in_models_library, flag_crs_in_models_with_resampling_library
from ._fileio import save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_models,
    save_intermediate_results,
    good_bits,
    maskpt,
    snr1,
    snr2,
    scale1,
    scale2,
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
    """
    Flag outliers in imaging data.

    input_models is expected to be a ModelLibrary

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelLibrary):
        on_disk = not in_memory
        input_models = ModelLibrary(input_models, on_disk=on_disk)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models

    if resample_data:
        # Start by creating resampled/mosaic images for
        # each group of exposures
        with input_models:
            example_model = input_models.borrow(0)
            output_path = make_output_path(basepath=example_model.meta.filename,
                        suffix='i2d')
            input_models.shelve(example_model, modify=False)
            del example_model
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
        median_wcs = resamp.output_wcs
        drizzled_models = resamp.do_drizzle(input_models)
    else:
        # for non-dithered data, the resampled image is just the original image
        drizzled_models = input_models
        with input_models:
            for i, model in enumerate(input_models):
                model.wht = build_driz_weight(
                    model,
                    weight_type=weight_type,
                    good_bits=good_bits)
                # copy for when saving median and input is a filename?
                if i == 0:
                    median_wcs = copy.deepcopy(model.meta.wcs)

    # Perform median combination on set of drizzled mosaics
    on_disk = not in_memory
    median_data = create_median_library(drizzled_models, maskpt, on_disk=on_disk)

    if save_intermediate_results:
        # make a median model
        with drizzled_models:
            example_model = drizzled_models.borrow(0)
            drizzled_models.shelve(example_model, modify=False)
            with datamodels.open(example_model) as dm0:
                median_model = datamodels.ImageModel(median_data)
                median_model.update(dm0)
                median_model.meta.wcs = median_wcs
            del example_model

        save_median(median_model, make_output_path, asn_id)
        del median_model

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    if resample_data:
        flag_crs_in_models_with_resampling_library(
            input_models,
            median_data,
            median_wcs,
            snr1,
            snr2,
            scale1,
            scale2,
            backg,
        )
    else:
        flag_crs_in_models_library(input_models, median_data, snr1)
    return input_models
