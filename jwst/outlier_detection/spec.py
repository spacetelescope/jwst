"""Perform outlier detection on spectra."""

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.utilities import record_step_status

from jwst.resample import resample_spec
from .utils import (
    flag_crs_in_models,
    flag_crs_in_models_with_resampling,
    median_with_resampling,
    median_without_resampling,
)

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
    make_output_path,
):
    """
    Flag outliers in slit-like spectroscopic data.

    Parameters
    ----------
    input_models : ModelContainer
        A container of data models.
    save_intermediate_results : bool
        If True, save intermediate results.
    good_bits : int
        Bit values indicating good pixels.
    maskpt : float
        The percentage of the mean weight to use as a threshold for masking.
    snr1 : float
        The signal-to-noise ratio threshold for first pass flagging, prior to smoothing.
    snr2 : float
        The signal-to-noise ratio threshold for secondary flagging, after smoothing.
    scale1 : float
        Scale factor used to scale the absolute derivative of the blot model for the first pass.
    scale2 : float
        Scale factor used to scale the absolute dervative of the blot model for the second pass.
    backg : float
        Scalar background level to add to the blotted image.
        Ignored if `input_model.meta.background.level` is not None but
        `input_model.meta.background.subtracted` is False.
    resample_data : bool
        If True, resample the data before detecting outliers.
    weight_type : str
        The type of weighting kernel to use when resampling.
        Options are 'ivm' or 'exptime'.
    pixfrac : float
        The pixel shrinkage factor to pass to drizzle.
    kernel : str
        The flux distribution kernel function to use when resampling.
    fillval : str
        The value to use in the output for pixels with no weight or flux
    make_output_path : function
        The functools.partial instance to pass to save_blot. Must be
        specified if save_blot is True.

    Returns
    -------
    ModelContainer
        The input models with outliers flagged.
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
        resamp = resample_spec.ResampleSpec(
            input_models,
            blendheaders=False,
            weight_type=weight_type,
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            good_bits=good_bits,
            enable_ctx=False,
            enable_var=False,
            compute_err="driz_err",
        )

        median_data, median_wcs, median_err = median_with_resampling(
            library,
            resamp,
            maskpt,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
            return_error=True,
        )
    else:
        median_data, median_wcs, median_err = median_without_resampling(
            library,
            maskpt,
            weight_type,
            good_bits,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
            return_error=True,
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
            make_output_path=make_output_path,
        )
    else:
        flag_crs_in_models(input_models, median_data, snr1, median_err=median_err)
    return input_models
