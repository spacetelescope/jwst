"""Utilities for outlier detection methods."""

import copy
from functools import partial
import numpy as np

from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample import compute_image_pixel_area
from stcal.resample.utils import build_driz_weight
from stcal.outlier_detection.utils import (
    compute_weight_threshold,
    gwcs_blot,
    flag_crs,
    flag_resampled_crs,
)
from stcal.outlier_detection.median import MedianComputer, nanmedian3D
from stdatamodels.jwst import datamodels
from . import _fileio

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel["DO_NOT_USE"]
OUTLIER = datamodels.dqflags.pixel["OUTLIER"]


def create_cube_median(cube_model, maskpt):
    """
    Compute the median over a cube of data.

    Parameters
    ----------
    cube_model : ~jwst.datamodels.CubeModel
        The input cube model.
    maskpt : float
        The percent threshold for masking bad data.

    Returns
    -------
    np.ndarray
        The median over the zeroth axis of the input cube.
    """
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)
    masked_cube = np.ma.masked_array(
        cube_model.data, np.less(cube_model.wht, weight_threshold)
    ).filled(np.nan)

    # not safe to use overwrite_input=True here because we are operating on model.data directly
    return nanmedian3D(masked_cube, overwrite_input=False)


def median_without_resampling(
    input_models,
    maskpt,
    weight_type,
    good_bits,
    save_intermediate_results=False,
    make_output_path=None,
    buffer_size=None,
    return_error=False,
):
    """
    Compute a median image without resampling.

    The median is performed across input exposures, for both
    imaging and spectral modes.

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.
    maskpt : float
        The weight threshold for masking out low weight pixels.
    weight_type : str
        The type of weighting to use when combining images. Options are:
        'ivm' (inverse variance) or 'exptime' (exposure time).
    good_bits : int
        The bit values that are considered good when determining the
        data quality of the input.
    save_intermediate_results : bool
        If True, save the drizzled models and median model to fits.
    make_output_path : function
        The functools.partial instance to pass to save_median. Must be
        specified if save_intermediate_results is True. Default None.
    buffer_size : int
        The size of chunk in bytes that will be read into memory when
        computing the median. This parameter has no effect if the input
        library has its on_disk attribute set to False.
    return_error : bool, optional
        If True, an approximate median error is computed alongside the
        median science image.

    Returns
    -------
    median_data : np.ndarray
        The median data array.
    median_wcs : gwcs.WCS
        A WCS corresponding to the median data.
    median_error : np.ndarray, optional
        A median error estimate, returned only if `return_error` is True.
    """
    in_memory = not input_models.on_disk
    ngroups = len(input_models)

    if save_intermediate_results:
        # create an empty image model for the median data
        median_model = datamodels.ImageModel(None)

    with input_models:
        for i in range(len(input_models)):
            drizzled_model = input_models.borrow(i)
            drizzled_data = drizzled_model.data.copy()
            if return_error:
                drizzled_err = drizzled_model.err.copy()
            else:
                drizzled_err = None
            weight = build_driz_weight(
                drizzled_model,
                weight_type=weight_type,
                good_bits=good_bits,
                flag_name_map=datamodels.dqflags.pixel,
            )
            if i == 0:
                median_wcs = copy.deepcopy(drizzled_model.meta.wcs)
                input_shape = (ngroups,) + drizzled_data.shape
                dtype = drizzled_data.dtype
                computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
                if return_error:
                    err_computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
                else:
                    err_computer = None
                if save_intermediate_results:
                    # update median model's meta with meta from the first model:
                    median_model.update(drizzled_model)
                    median_model.meta.wcs = median_wcs

            weight_threshold = compute_weight_threshold(weight, maskpt)
            drizzled_data[weight < weight_threshold] = np.nan
            computer.append(drizzled_data, i)
            if return_error:
                drizzled_err[weight < weight_threshold] = np.nan
                err_computer.append(drizzled_err, i)

            input_models.shelve(drizzled_model, i, modify=False)
            del drizzled_model

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()
    if return_error:
        median_err = err_computer.evaluate()
    else:
        median_err = None

    if save_intermediate_results:
        # Save median model to fits
        median_model.data = median_data
        if return_error:
            median_model.err = median_err
        _fileio.save_median(median_model, make_output_path)

    if return_error:
        return median_data, median_wcs, median_err
    else:
        return median_data, median_wcs


def median_with_resampling(
    input_models,
    resamp,
    maskpt,
    save_intermediate_results=False,
    make_output_path=None,
    buffer_size=None,
    return_error=False,
):
    """
    Compute a median image with resampling.

    The median is performed across resampled groups, for both imaging
    and spectral modes.

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.
    resamp : resample.resample.ResampleImage object
        The controlling object for the resampling process.
    maskpt : float
        The weight threshold for masking out low weight pixels.
    save_intermediate_results : bool
        If True, save the drizzled models and median model to fits.
    make_output_path : function
        The functools.partial instance to pass to save_median. Must be
        specified if save_intermediate_results is True. Default None.
    buffer_size : int
        The size of chunk in bytes that will be read into memory when
        computing the median. This parameter has no effect if the input
        library has its on_disk attribute set to False.
    return_error : bool, optional
        If True, an approximate median error is computed alongside the
        median science image.

    Returns
    -------
    median_data : np.ndarray
        The median data array.
    median_wcs : gwcs.WCS
        A WCS corresponding to the median data.
    median_error : np.ndarray, None, optional
        A median error estimate, returned only if `return_error` is `True`.
        If ``resamp.compute_err`` is not set to "driz_err", `None` will be
        returned.
    """
    in_memory = not input_models.on_disk
    indices_by_group = list(input_models.group_indices.values())
    ngroups = len(indices_by_group)
    median_err = None

    eval_med_err = False
    if return_error:
        if resamp.compute_err == "driz_err":
            eval_med_err = True
        else:
            log.warning(
                "Returning median_error has been disabled since input "
                "'resamp' object does not have 'compute_err' attribute set to "
                "'driz_err'."
            )

    if save_intermediate_results:
        # create an empty image model for the median data
        median_model = datamodels.ImageModel(None)

    for i, indices in enumerate(indices_by_group):
        drizzled_model = resamp.resample_group(indices)

        if save_intermediate_results:
            # write the drizzled model to file
            _fileio.save_drizzled(drizzled_model, make_output_path)

        if i == 0:
            median_wcs = resamp.output_wcs
            input_shape = (ngroups,) + drizzled_model.data.shape
            dtype = drizzled_model.data.dtype
            computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
            if eval_med_err:
                err_computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
            else:
                err_computer = None
            if save_intermediate_results:
                # update median model's meta with meta from the first model:
                median_model.update(drizzled_model)
                median_model.meta.wcs = median_wcs

        weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
        drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan
        computer.append(drizzled_model.data, i)
        if eval_med_err:
            drizzled_model.err[drizzled_model.wht < weight_threshold] = np.nan
            err_computer.append(drizzled_model.err, i)
        del drizzled_model

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()
    if eval_med_err:
        median_err = err_computer.evaluate()

    if save_intermediate_results:
        # Save median model to fits
        median_model.data = median_data
        if eval_med_err:
            median_model.err = median_err
        # drizzled model already contains asn_id
        make_output_path = partial(make_output_path, asn_id=None)
        _fileio.save_median(median_model, make_output_path)

    if return_error:
        return median_data, median_wcs, median_err
    else:
        return median_data, median_wcs


def flag_crs_in_models(input_models, median_data, snr1, median_err=None):
    """
    Flag outliers in all input models without resampling.

    Parameters
    ----------
    input_models : ModelContainer
        The input datamodels.
    median_data : np.ndarray
        The median data array.
    snr1 : float
        The signal-to-noise ratio threshold for flagging outliers.
    median_err : np.ndarray, optional
        The error array corresponding to the median data. If not provided,
        the error array stored the input model `err` extension will be used.
    """
    for image in input_models:
        # dq flags will be updated in-place
        flag_model_crs(image, median_data, snr1, median_err=median_err)


def flag_resampled_model_crs(
    input_model,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    median_err=None,
    save_blot=False,
    make_output_path=None,
):
    """
    Flag outliers in a resampled model, updating DQ array in place.

    Parameters
    ----------
    input_model : ~jwst.datamodels.DataModel
        The input datamodel.
    median_data : np.ndarray
        The median data array.
    median_wcs : gwcs.WCS
        A WCS corresponding to the median data.
    snr1 : float
        The signal-to-noise ratio threshold for first pass flagging, prior to smoothing.
    snr2 : float
        The signal-to-noise ratio threshold for secondary flagging, after smoothing.
    scale1 : float
        Scale factor used to scale the absolute derivative of the blot model for the first pass.
    scale2 : float
        Scale factor used to scale the absolute derivative of the blot model for the second pass.
    backg : float
        Scalar background level to add to the blotted image.
        Ignored if `input_model.meta.background.level` is not None but
        `input_model.meta.background.subtracted` is False.
    median_err : np.ndarray, optional
        The error array corresponding to the median data. If not provided,
        the error array stored the input model `err` extension will be used.
    save_blot : bool
        If True, save the blotted image to fits.
    make_output_path : function
        The functools.partial instance to pass to save_blot. Must be
        specified if save_blot is True.
    """
    if "SPECTRAL" not in input_model.meta.wcs.output_frame.axes_type:
        input_pixflux_area = input_model.meta.photometry.pixelarea_steradians
        # Set array shape, needed to compute image pixel area
        input_model.meta.wcs.array_shape = input_model.shape
        input_pixel_area = compute_image_pixel_area(input_model.meta.wcs)
        pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
    else:
        pix_ratio = 1.0

    blot = gwcs_blot(
        median_data,
        median_wcs,
        input_model.data.shape,
        input_model.meta.wcs,
        pix_ratio,
        fillval=np.nan,
    )
    if median_err is not None:
        blot_err = gwcs_blot(
            median_err,
            median_wcs,
            input_model.data.shape,
            input_model.meta.wcs,
            pix_ratio,
            fillval=np.nan,
        )
    else:
        blot_err = None
    if save_blot:
        _fileio.save_blot(input_model, blot, blot_err, make_output_path)

    # dq flags will be updated in-place
    _flag_resampled_model_crs(input_model, blot, blot_err, snr1, snr2, scale1, scale2, backg)


def _flag_resampled_model_crs(
    input_model,
    blot,
    blot_err,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    """
    Flag outliers via comparison to a blotted image and update the DQ array in place.

    Parameters
    ----------
    input_model : ~jwst.datamodels.DataModel
        The input datamodel.
    blot : np.ndarray
        The blotted data array.
    blot_err : np.ndarray
        The blotted error array.
    snr1 : float
        The signal-to-noise ratio threshold for first pass flagging, prior to smoothing.
    snr2 : float
        The signal-to-noise ratio threshold for secondary flagging, after smoothing.
    scale1 : float
        Scale factor used to scale the absolute derivative of the blot model for the first pass.
    scale2 : float
        Scale factor used to scale the absolute derivative of the blot model for the second pass.
    backg : float
        Scalar background level to add to the blotted image.
        Ignored if `input_model.meta.background.level` is not None but
        `input_model.meta.background.subtracted` is False.
    """
    if (
        input_model.meta.background.subtracted is False
        and input_model.meta.background.level is not None
    ):
        backg = input_model.meta.background.level
        log.debug(f"Adding background level {backg} to blotted image")

    if blot_err is not None:
        err_to_use = blot_err
    else:
        err_to_use = input_model.err
    cr_mask = flag_resampled_crs(
        input_model.data, err_to_use, blot, snr1, snr2, scale1, scale2, backg
    )

    # update the dq flags in-place
    input_model.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)
    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")

    # Make sure all data, error, and variance arrays have
    # matching NaNs and DQ flags
    match_nans_and_flags(input_model)


def flag_crs_in_models_with_resampling(
    input_models,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    median_err=None,
    save_blot=False,
    make_output_path=None,
):
    """
    Flag outliers in all input models, with resampling, modifying DQ array in place.

    Parameters
    ----------
    input_models : ~jwst.datamodels.ModelContainer
        The input datamodels.
    median_data : np.ndarray
        The median data array.
    median_wcs : gwcs.WCS
        A WCS corresponding to the median data.
    snr1 : float
        The signal-to-noise ratio threshold for first pass flagging, prior to smoothing.
    snr2 : float
        The signal-to-noise ratio threshold for secondary flagging, after smoothing.
    scale1 : float
        Scale factor used to scale the absolute derivative of the blot model for the first pass.
    scale2 : float
        Scale factor used to scale the absolute derivative of the blot model for the second pass.
    backg : float
        Scalar background level to add to the blotted image.
        Ignored if `input_model.meta.background.level` is not None but
        `input_model.meta.background.subtracted` is False.
    median_err : np.ndarray, optional
        The error array corresponding to the median data. If not provided,
        the error array stored the input model `err` extension will be used.
    save_blot : bool
        If True, save the blotted image to fits.
    make_output_path : function
        The functools.partial instance to pass to save_blot. Must be
        specified if save_blot is True.
    """
    for image in input_models:
        flag_resampled_model_crs(
            image,
            median_data,
            median_wcs,
            snr1,
            snr2,
            scale1,
            scale2,
            backg,
            median_err=median_err,
            save_blot=save_blot,
            make_output_path=make_output_path,
        )


def flag_model_crs(image, blot, snr, median_err=None):
    """
    Flag outliers in a model.

    Parameters
    ----------
    image : ~jwst.datamodels.DataModel
        The input datamodel.
    blot : np.ndarray
        The blotted data array.
    snr : float
        The signal-to-noise ratio threshold for flagging outliers.
    median_err : np.ndarray, optional
        The error array corresponding to the median data. If not provided,
        the error array stored the input model `err` extension will be used.
    """
    if median_err is not None:
        error_to_use = median_err
    else:
        error_to_use = image.err
    cr_mask = flag_crs(image.data, error_to_use, blot, snr)

    # Update dq array in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    # Make sure all data, error, and variance arrays have
    # matching NaNs and DQ flags
    match_nans_and_flags(image)

    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")
