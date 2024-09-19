"""
Submodule for performing outlier detection on imaging data.
"""

import logging
import os

from jwst.datamodels import ModelLibrary
from jwst.resample import resample
from jwst.stpipe.utilities import record_step_status

from .utils import flag_model_crs, flag_resampled_model_crs, drizzle_and_median

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
        input_models = ModelLibrary(input_models, on_disk=not in_memory)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models
    
    with input_models:
        example_model = input_models.borrow(0)
        output_path = make_output_path(basepath=example_model.meta.filename,
                    suffix='')
        input_models.shelve(example_model, modify=False)
        del example_model
    output_path = os.path.dirname(output_path)
        
    # Start by creating resampled/mosaic images for
    # each group of exposures
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

    median_data, median_wcs = drizzle_and_median(input_models,
                                                 resamp,
                                                 maskpt,
                                                 resample_data=resample_data,
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
