"""
Class definition for performing outlier detection on spectra.

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

from stdatamodels.jwst import datamodels

from ..resample import resample_spec, resample_utils
from .utils import _convert_inputs, _detect_outliers, _remove_file, create_median

import logging
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
    in_memory,
    asn_id,
    make_output_path,
):
    """Flag outlier pixels in DQ of input images."""
    input_models = _convert_inputs(input_models, good_bits, weight_type)

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
            in_memory=in_memory,
            asn_id=asn_id,
        )
        drizzled_models = resamp.do_drizzle(input_models)
        if save_intermediate_results:
            for model in drizzled_models:
                model.meta.filename = make_output_path(
                    basepath=model.meta.filename,
                    suffix="_outlier_s2d.fits",
                )
                log.info("Writing out resampled spectra...")
                model.save(model.meta.filename)
    else:
        drizzled_models = input_models
        for i in range(len(input_models)):
            drizzled_models[i].wht = resample_utils.build_driz_weight(
                input_models[i],
                weight_type=weight_type,
                good_bits=good_bits)

    # Initialize intermediate products used in the outlier detection
    median_model = datamodels.ImageModel(drizzled_models[0].data.shape)
    median_model.meta = drizzled_models[0].meta
    median_model.meta.filename = make_output_path(
        basepath=input_models[0].meta.filename,
        suffix='median'
    )

    # Perform median combination on set of drizzled mosaics
    # create_median should be called as a method from parent class
    median_model.data = create_median(drizzled_models, maskpt)

    if save_intermediate_results:
        log.info("Writing out MEDIAN image to: {}".format(
                 median_model.meta.filename))
        median_model.save(median_model.meta.filename)
    else:
        # since we're not saving intermediate results if the drizzled models
        # were written to disk, remove them
        if not in_memory:
            for fn in drizzled_models._models:
                _remove_file(fn)
                log.info(f"Removing file {fn}")

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

    del median_model
