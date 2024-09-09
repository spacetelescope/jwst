"""
Submodule for performing outlier detection on spectra.
"""
import copy
import os

from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.utilities import record_step_status

from ..resample import resample_spec, resample_utils
from .utils import create_median, flag_crs_in_models, flag_crs_in_models_with_resampling
from ._fileio import remove_file

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
    asn_id,
    make_output_path,
):
    """
    Flag outliers in spec data.

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelContainer):
        input_models = ModelContainer(input_models, save_open=in_memory)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models

    if resample_data is True:
        # Start by creating resampled/mosaic images for
        #  each group of exposures
        output_path = make_output_path(
            basepath=input_models[0].meta.filename, suffix='')
        output_path = os.path.dirname(output_path)
        resamp = resample_spec.ResampleSpecData(
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
        )
        median_wcs = resamp.output_wcs

        # convert to library for resample
        # for compatibility with image3 pipeline
        library = ModelLibrary(input_models, on_disk=False)
        drizzled_models = resamp.do_drizzle(library)

    else:
        drizzled_models = ModelLibrary(input_models)
        with drizzled_models:
            for i, model in enumerate(drizzled_models):
                model.wht = resample_utils.build_driz_weight(
                    input_models[i],
                    weight_type=weight_type,
                    good_bits=good_bits)
                drizzled_models.shelve(model)
        # copy for when saving median and input is a filename?
        median_wcs = copy.deepcopy(input_models[0].meta.wcs)

    # Perform median combination on set of drizzled mosaics
    # create_median should be called as a method from parent class
    median_data = create_median(drizzled_models, maskpt)

    if save_intermediate_results:
        # Initialize intermediate products used in the outlier detection
        median_model = datamodels.ImageModel(median_data)
        with drizzled_models:
            example_model = drizzled_models.borrow(0)
            drizzled_models.shelve(example_model, 0, modify=False)
        median_model.meta = example_model.meta
        median_model.meta.filename = make_output_path(
            basepath=input_models[0].meta.filename,
            suffix='median'
        )

        log.info("Writing out MEDIAN image to: {}".format(
                 median_model.meta.filename))
        median_model.save(median_model.meta.filename)
        del median_model
    else:
        # since we're not saving intermediate results if the drizzled models
        # were written to disk, remove them
        if not in_memory:
            for fn in drizzled_models._members:
                remove_file(fn["expname"])

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
            save_blot=save_intermediate_results,
            make_output_path=make_output_path
        )
    else:
        flag_crs_in_models(input_models,
                           median_data,
                           snr1)
    return input_models
