"""
The ever-present utils sub-module. A home for all...
"""
import warnings

import numpy as np

from jwst.resample.resample import compute_image_pixel_area
from stcal.outlier_detection.utils import compute_weight_threshold, gwcs_blot, flag_crs, flag_resampled_crs
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    median = np.nanmedian(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        axis=0,
    )
    return median


def create_median(resampled_models, maskpt):
    """Create a median image from the singly resampled images.
    """
    log.info("Computing median")

    weight_thresholds = compute_weight_threshold_container(resampled_models, maskpt)

    # Now, set up buffered access to all input models
    resampled_models.set_buffer(1.0)  # Set buffer at 1Mb
    resampled_sections = resampled_models.get_sections()
    median_image = np.empty((resampled_models.imrows, resampled_models.imcols),
                            resampled_models.imtype)
    median_image[:] = np.nan  # initialize with NaNs

    for (resampled_sci, resampled_weight, (row1, row2)) in resampled_sections:
        # Create a mask for each input image, masking out areas where there is
        # no data or the data has very low weight
        badmasks = []
        for weight, weight_threshold in zip(resampled_weight, weight_thresholds):
            badmask = np.less(weight, weight_threshold)
            log.debug("Percentage of pixels with low weight: {}".format(
                np.sum(badmask) / len(weight.flat) * 100))
            badmasks.append(badmask)

        # Fill resampled_sci array with nan's where mask values are True
        for f1, f2 in zip(resampled_sci, badmasks):
            for elem1, elem2 in zip(f1, f2):
                elem1[elem2] = np.nan

        del badmasks

        # For a stack of images with "bad" data replaced with Nan
        # use np.nanmedian to compute the median.
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore",
                                    message="All-NaN slice encountered",
                                    category=RuntimeWarning)
            median_image[row1:row2] = np.nanmedian(resampled_sci, axis=0)
        del resampled_sci, resampled_weight

    return median_image


def compute_weight_threshold_container(resampled_models, maskpt):
    '''
    Compute weight means without keeping datamodels for each input open

    Parameters
    ----------
    resampled_models : ~jwst.datamodels.ModelContainer
        The input data models.

    maskpt : float
        The percentage of the mean weight to use as a threshold for masking.

    Returns
    -------
    list
        The weight thresholds for each integration.
    '''

    # Start by ensuring that the ModelContainer does NOT open and keep each datamodel
    ropen_orig = resampled_models._return_open
    resampled_models._return_open = True
    # keep track of resulting computation for each input resampled datamodel
    weight_thresholds = []
    # For each model, compute the bad-pixel threshold from the weight arrays
    for resampled in resampled_models:
        weight = resampled.wht
        weight_threshold = compute_weight_threshold(weight, maskpt)
        weight_thresholds.append(weight_threshold)
        # close and delete the model, just to explicitly try to keep the memory as clean as possible
        resampled.close()
        del resampled
    # Reset ModelContainer attribute to original value
    resampled_models._return_open = ropen_orig
    return weight_thresholds


def flag_crs_in_models(
    input_models,
    median_data,
    snr1,
):
    for image in input_models:
        # dq flags will be updated in-place
        flag_model_crs(image, median_data, snr1)


def flag_crs_in_models_with_resampling(
    input_models,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    for image in input_models:
        if 'SPECTRAL' not in image.meta.wcs.output_frame.axes_type:
            input_pixflux_area = image.meta.photometry.pixelarea_steradians
            # Set array shape, needed to compute image pixel area
            image.meta.wcs.array_shape = image.shape
            input_pixel_area = compute_image_pixel_area(image.meta.wcs)
            pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
        else:
            pix_ratio = 1.0

        blot = gwcs_blot(median_data, median_wcs, image.data.shape, image.meta.wcs, pix_ratio)
        # dq flags will be updated in-place
        flag_resampled_model_crs(image, blot, snr1, snr2, scale1, scale2, backg)


def flag_resampled_model_crs(
    image,
    blot,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    """
    Flag crs in image based on a resampled (and blotted) data (blot).
    """
    # If the datamodel has a measured background that has not been subtracted
    # use it instead of the user provided backg.
    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (image.meta.background.subtracted is False and
            image.meta.background.level is not None):
        backg = image.meta.background.level
        log.debug(f"Adding background level {backg} to blotted image")

    cr_mask = flag_resampled_crs(image.data, image.err, blot, snr1, snr2, scale1, scale2, backg)

    # update the dq flags in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")


def flag_model_crs(image, blot, snr):
    cr_mask = flag_crs(image.data, image.err, blot, snr)
    # update dq array in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)
    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")
