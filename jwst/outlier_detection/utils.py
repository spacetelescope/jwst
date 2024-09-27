"""
The ever-present utils sub-module. A home for all...
"""

import copy
from functools import partial
import numpy as np

from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample import compute_image_pixel_area
from jwst.resample.resample_utils import build_driz_weight
from stcal.outlier_detection.utils import compute_weight_threshold, gwcs_blot, flag_crs, flag_resampled_crs
from stcal.outlier_detection.median import MedianComputer, nanmedian3D
from stdatamodels.jwst import datamodels
from . import _fileio

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    # not safe to use overwrite_input=True here because we are operating on model.data directly
    return nanmedian3D(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        overwrite_input=False)


def median_without_resampling(input_models,
                              maskpt,
                              weight_type,
                              good_bits,
                              save_intermediate_results=False,
                              make_output_path=None,
                              buffer_size=None):
    """
    Shared code between imaging and spec modes for resampling and median computation

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
        if True, save the drizzled models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median. Must be 
        specified if save_intermediate_results is True. Default None.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False.
    """
    in_memory = not input_models._on_disk
    ngroups = len(input_models)

    with input_models:
        for i in range(len(input_models)):

            drizzled_model = input_models.borrow(i)
            drizzled_model.wht = build_driz_weight(drizzled_model,
                                                    weight_type=weight_type,
                                                    good_bits=good_bits)
            median_wcs = copy.deepcopy(drizzled_model.meta.wcs)
            input_models.shelve(drizzled_model, i, modify=True)

            if save_intermediate_results:
                # write the drizzled model to file
                _fileio.save_drizzled(drizzled_model, make_output_path)

            if i == 0:
                input_shape = (ngroups,)+drizzled_model.data.shape
                dtype = drizzled_model.data.dtype
                computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)

            weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
            drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan
            computer.append(drizzled_model.data, i)

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()

    if save_intermediate_results:
        # Save median model to fits
        median_model = datamodels.ImageModel(median_data)
        median_model.update(drizzled_model)
        median_model.meta.wcs = median_wcs
        _fileio.save_median(median_model, make_output_path)
    del drizzled_model

    return median_data, median_wcs


def median_with_resampling(input_models,
                           resamp,
                           maskpt,
                           save_intermediate_results=False,
                           make_output_path=None,
                           buffer_size=None):
    """
    Shared code between imaging and spec modes for resampling and median computation

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.

    resamp : resample.resample.ResampleData object
        The controlling object for the resampling process.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    save_intermediate_results : bool
        if True, save the drizzled models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median. Must be 
        specified if save_intermediate_results is True. Default None.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False.
    """
    if not resamp.single:
        raise ValueError("median_with_resampling should only be used for resample_many_to_many")
    
    in_memory = not input_models._on_disk
    indices_by_group = list(input_models.group_indices.values())
    ngroups = len(indices_by_group)

    with input_models:
        for i, indices in enumerate(indices_by_group):

            median_wcs = resamp.output_wcs
            drizzled_model = resamp.resample_group(input_models, indices)

            if save_intermediate_results:
                # write the drizzled model to file
                _fileio.save_drizzled(drizzled_model, make_output_path)

            if i == 0:
                input_shape = (ngroups,)+drizzled_model.data.shape
                dtype = drizzled_model.data.dtype
                computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)

            weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
            drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan
            computer.append(drizzled_model.data, i)

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()

    if save_intermediate_results:
        # Save median model to fits
        median_model = datamodels.ImageModel(median_data)
        median_model.update(drizzled_model)
        median_model.meta.wcs = median_wcs
        # drizzled model already contains asn_id
        make_output_path = partial(make_output_path, asn_id=None)
        _fileio.save_median(median_model, make_output_path)
    del drizzled_model

    return median_data, median_wcs


def flag_crs_in_models(
    input_models,
    median_data,
    snr1,
):
    for image in input_models:
        # dq flags will be updated in-place
        flag_model_crs(image, median_data, snr1)
    

def flag_resampled_model_crs(
    input_model,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    save_blot=False,
    make_output_path=None,
):
    if 'SPECTRAL' not in input_model.meta.wcs.output_frame.axes_type:
        input_pixflux_area = input_model.meta.photometry.pixelarea_steradians
        # Set array shape, needed to compute image pixel area
        input_model.meta.wcs.array_shape = input_model.shape
        input_pixel_area = compute_image_pixel_area(input_model.meta.wcs)
        pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
    else:
        pix_ratio = 1.0

    blot = gwcs_blot(median_data, median_wcs, input_model.data.shape, input_model.meta.wcs, pix_ratio)
    if save_blot:
        _fileio.save_blot(input_model, blot, make_output_path)
    # dq flags will be updated in-place
    _flag_resampled_model_crs(input_model, blot, snr1, snr2, scale1, scale2, backg)


def _flag_resampled_model_crs(
    input_model,
    blot,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    # If the datamodel has a measured background that has not been subtracted
    # use it instead of the user provided backg.
    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (input_model.meta.background.subtracted is False and
            input_model.meta.background.level is not None):
        backg = input_model.meta.background.level
        log.debug(f"Adding background level {backg} to blotted image")

    cr_mask = flag_resampled_crs(input_model.data, input_model.err, blot, snr1, snr2, scale1, scale2, backg)

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
    save_blot=False,
    make_output_path=None,
):
    for image in input_models:
        flag_resampled_model_crs(image,
                                 median_data,
                                 median_wcs,
                                 snr1,
                                 snr2,
                                 scale1,
                                 scale2,
                                 backg,
                                 save_blot=save_blot,
                                 make_output_path=make_output_path)


def flag_model_crs(image, blot, snr):
    cr_mask = flag_crs(image.data, image.err, blot, snr)

    # Update dq array in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    # Make sure all data, error, and variance arrays have
    # matching NaNs and DQ flags
    match_nans_and_flags(image)

    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")
