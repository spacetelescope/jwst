"""Perform outlier detection on spectra."""

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.utilities import record_step_status

from ..resample import resample_spec
from .utils import (flag_crs_in_models,
                    flag_crs_in_models_with_resampling,
                    median_with_resampling,
                    median_without_resampling)

import logging
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
    """Flag outliers in spec data.

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelContainer):
        input_models = ModelContainer(input_models)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models

    # convert to library for resample
    # for compatibility with image3 pipeline
    library = ModelLibrary(input_models, on_disk=False)

    if resample_data is True:
        # Start by creating resampled/mosaic images for
        #  each group of exposures
        resamp = resample_spec.ResampleSpecData(
            input_models,
            single=True,
            blendheaders=False,
            wht_type=weight_type,
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            good_bits=good_bits,
        )

        median_data, median_wcs, median_err = median_with_resampling(
            library,
            resamp,
            maskpt,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
            return_error=True)
    else:
        median_data, median_wcs, median_err = median_without_resampling(
            library,
            maskpt,
            weight_type,
            good_bits,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
            return_error=True
        )

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    if resample_data:
        flag_crs_in_models_with_resampling(
            input_models,
            median_data,
            median_wcs,
            snr1,
            snr2,
            scale1,
            scale2,
            backg,
            median_err=median_err,
            save_blot=save_intermediate_results,
            make_output_path=make_output_path
        )
    else:
        flag_crs_in_models(input_models,
                           median_data,
                           snr1,
                           median_err=median_err)
    return input_models
