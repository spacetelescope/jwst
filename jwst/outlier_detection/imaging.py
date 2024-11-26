"""
Submodule for performing outlier detection on imaging data.
"""

import logging

from jwst.datamodels import ModelLibrary
from jwst.resample import resample
from jwst.stpipe.utilities import record_step_status

from .utils import (flag_model_crs,
                    flag_resampled_model_crs,
                    median_without_resampling,
                    median_with_resampling)

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
    in_memory,
    make_output_path,
):
    """
    Flag outliers in imaging data.

    input_models is expected to be a ModelLibrary

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelLibrary):
        input_models = ModelLibrary(input_models, on_disk=not in_memory)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models
        
    if resample_data:
        resamp = resample.ResampleData(
            input_models,
            single=True,
            blendheaders=False,
            wht_type=weight_type,
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            good_bits=good_bits,
        )
        median_data, median_wcs = median_with_resampling(input_models,
                                                    resamp,
                                                    maskpt,
                                                    save_intermediate_results=save_intermediate_results,
                                                    make_output_path=make_output_path,)
    else:
        median_data, median_wcs = median_without_resampling(input_models,
                                                    maskpt,
                                                    weight_type,
                                                    good_bits,
                                                    save_intermediate_results=save_intermediate_results,
                                                    make_output_path=make_output_path,)


    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    with input_models:
        for image in input_models:
            if resample_data:
                flag_resampled_model_crs(image,
                                         median_data,
                                         median_wcs,
                                         snr1,
                                         snr2,
                                         scale1,
                                         scale2,
                                         backg,
                                         save_blot=save_intermediate_results,
                                         make_output_path=make_output_path)
            else:
                flag_model_crs(image, median_data, snr1)
            input_models.shelve(image, modify=True)

    return input_models
