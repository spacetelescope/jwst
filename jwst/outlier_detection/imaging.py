"""
Submodule for performing outlier detection on imaging data.
"""

import copy
import logging
import os

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight
from jwst.stpipe.utilities import record_step_status

from .utils import create_median, flag_crs_in_models, flag_crs_in_models_with_resampling
from ._fileio import remove_file, save_median

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

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelContainer):
        input_models = ModelContainer(input_models, save_open=in_memory)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models

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
        median_wcs = resamp.output_wcs
        drizzled_models = resamp.do_drizzle(input_models)
    else:
        # for non-dithered data, the resampled image is just the original image
        drizzled_models = input_models
        for i in range(len(input_models)):
            drizzled_models[i].wht = build_driz_weight(
                input_models[i],
                weight_type=weight_type,
                good_bits=good_bits)
        # copy for when saving median and input is a filename?
        median_wcs = copy.deepcopy(input_models[0].meta.wcs)

    # Perform median combination on set of drizzled mosaics
    median_data = create_median(drizzled_models, maskpt)

    if save_intermediate_results:
        # make a median model
        with datamodels.open(drizzled_models[0]) as dm0:
            median_model = datamodels.ImageModel(median_data)
            median_model.update(dm0)
            median_model.meta.wcs = median_wcs

        save_median(median_model, make_output_path, asn_id)
        del median_model
    else:
        # since we're not saving intermediate results if the drizzled models
        # were written to disk, remove them
        if not in_memory:
            for fn in drizzled_models._models:
                remove_file(fn)

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
        )
    else:
        flag_crs_in_models(input_models, median_data, snr1)
    return input_models
