"""
The ever-present utils sub-module. A home for all...
"""
import os
import warnings

import numpy as np
from astropy.stats import sigma_clip
from drizzle.cdrizzle import tblot
from scipy import ndimage

from jwst.datamodels import ModelContainer
from jwst.resample.resample_utils import build_driz_weight, calc_gwcs_pixmap
from jwst.resample.resample import compute_image_pixel_area
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# flags used in flag_cr_update_model
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']


# TODO stcal
def compute_weight_threshold(weight, maskpt):
    '''
    Compute the weight threshold for a single image or cube.

    Parameters
    ----------
    weight : numpy.ndarray
        The weight array

    maskpt : float
        The percentage of the mean weight to use as a threshold for masking.

    Returns
    -------
    float
        The weight threshold for this integration.
    '''

    # necessary in order to assure that mask gets applied correctly
    if hasattr(weight, '_mask'):
        del weight._mask
    mask_zero_weight = np.equal(weight, 0.)
    mask_nans = np.isnan(weight)
    # Combine the masks
    weight_masked = np.ma.array(weight, mask=np.logical_or(
        mask_zero_weight, mask_nans))
    # Sigma-clip the unmasked data
    weight_masked = sigma_clip(weight_masked, sigma=3, maxiters=5)
    mean_weight = np.mean(weight_masked)
    # Mask pixels where weight falls below maskpt percent
    weight_threshold = mean_weight * maskpt
    return weight_threshold


# TODO stcal as a "median computer"
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


# TODO stcal utils
def _abs_deriv(array):
    """Take the absolute derivate of a numpy array."""
    tmp = np.zeros(array.shape, dtype=np.float64)
    out = np.zeros(array.shape, dtype=np.float64)

    tmp[1:, :] = array[:-1, :]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:-1, :] = array[1:, :]
    tmp, out = _absolute_subtract(array, tmp, out)

    tmp[:, 1:] = array[:, :-1]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:, :-1] = array[:, 1:]
    tmp, out = _absolute_subtract(array, tmp, out)

    return out


# TODO stcal utils
def _absolute_subtract(array, tmp, out):
    tmp = np.abs(array - tmp)
    out = np.maximum(tmp, out)
    tmp = tmp * 0.
    return tmp, out


# TODO stcal utils
def flag_cr(
    sci_data,
    sci_err,
    blot_data,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    resample_data,
):
    """Masks outliers in science image by updating DQ in-place

    Mask blemishes in dithered data by comparing a science image
    with a model image and the derivative of the model image.

    Parameters
    ----------

    FIXME: update these

    sci_image : ~jwst.datamodels.ImageModel
        the science data. Can also accept a CubeModel, but only if
        resample_data is False

    blot_array : np.ndarray
        the blotted median image of the dithered science frames.

    snr : str
        Signal-to-noise ratio

    scale : str
        scaling factor applied to the derivative

    backg : float
        Background value (scalar) to subtract

    resample_data : bool
        Boolean to indicate whether blot_image is created from resampled,
        dithered data or not

    Notes
    -----
    Accepting a CubeModel for sci_image and blot_image with resample_data=True
    appears to be a relatively simple extension, as the only thing that explicitly
    relies on the dimensionality is the kernel, which could be generalized.
    However, this is not currently needed, as CubeModels are only passed in for
    TSO data, where resampling is always False.
    """
    err_data = np.nan_to_num(sci_err)

    # create the outlier mask
    if resample_data:  # dithered outlier detection
        blot_deriv = _abs_deriv(blot_data)
        diff_noise = np.abs(sci_data - blot_data - backg)

        # Create a boolean mask based on a scaled version of
        # the derivative image (dealing with interpolating issues?)
        # and the standard n*sigma above the noise
        threshold1 = scale1 * blot_deriv + snr1 * err_data
        mask1 = np.greater(diff_noise, threshold1)

        # Smooth the boolean mask with a 3x3 boxcar kernel
        kernel = np.ones((3, 3), dtype=int)
        mask1_smoothed = ndimage.convolve(mask1, kernel, mode='nearest')

        # Create a 2nd boolean mask based on the 2nd set of
        # scale and threshold values
        threshold2 = scale2 * blot_deriv + snr2 * err_data
        mask2 = np.greater(diff_noise, threshold2)

        # Final boolean mask
        cr_mask = mask1_smoothed & mask2

    else:  # stack outlier detection
        diff_noise = np.abs(sci_data - blot_data)

        # straightforward detection of outliers for non-dithered data since
        # err_data includes all noise sources (photon, read, and flat for baseline)
        cr_mask = np.greater(diff_noise, snr1 * err_data)

    return cr_mask


# TODO stcal? (would require moving calc_gwcs_pixmap)
# FIXME (or fixed) interp and sinscl were "options" only when provided
# as part of the step spec (which becomes outlierpars). As neither was
# in the spec (and providing unknown arguments causes an error), these
# were never configurable and always defaulted to linear and 1.0
def gwcs_blot(median_data, median_wcs, blot_data, blot_wcs, pix_ratio):
    """
    Resample the output/resampled image to recreate an input image based on
    the input image's world coordinate system

    Parameters
    ----------
    median_model : `~stdatamodels.jwst.datamodels.JwstDataModel`

    blot_img : datamodel
        Datamodel containing header and WCS to define the 'blotted' image
    """
    # Compute the mapping between the input and output pixel coordinates
    pixmap = calc_gwcs_pixmap(blot_wcs, median_wcs, blot_data.shape)
    log.debug("Pixmap shape: {}".format(pixmap[:, :, 0].shape))
    log.debug("Sci shape: {}".format(blot_data.shape))
    log.info('Blotting {} <-- {}'.format(blot_data.shape, median_data.shape))

    outsci = np.zeros(blot_data.shape, dtype=np.float32)

    # Currently tblot cannot handle nans in the pixmap, so we need to give some
    # other value.  -1 is not optimal and may have side effects.  But this is
    # what we've been doing up until now, so more investigation is needed
    # before a change is made.  Preferably, fix tblot in drizzle.
    pixmap[np.isnan(pixmap)] = -1
    tblot(median_data, pixmap, outsci, scale=pix_ratio, kscale=1.0,
          interp='linear', exptime=1.0, misval=0.0, sinscl=1.0)

    return outsci


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


def _detect_outliers(
    input_models,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    resample_data,
):
    for image in input_models:
        if resample_data:
            if 'SPECTRAL' not in image.meta.wcs.output_frame.axes_type:
                input_pixflux_area = image.meta.photometry.pixelarea_steradians
                input_pixel_area = compute_image_pixel_area(image.meta.wcs)
                pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
            else:
                pix_ratio = 1.0

            blot = gwcs_blot(median_data, median_wcs, image.data, image.meta.wcs, pix_ratio)
        else:
            blot = median_data

        # dq flags will be updated in-place
        flag_cr_update_model(image, blot, snr1, snr2, scale1, scale2, backg, resample_data)


def flag_cr_update_model(
    image,
    blot,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    resample_data,
):
    # previous versions of the code generated all of the blot models and wrote
    # them to disk. Now the code generates the blot model(s) only when needed
    # and no longer needs to write them to disk. We could re-introduce saving
    # of blot models here.

    # If the datamodel has a measured background that has not been subtracted
    # use it instead of the user provided backg.
    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (image.meta.background.subtracted is False and
            image.meta.background.level is not None):
        backg = image.meta.background.level
        log.debug(f"Adding background level {backg} to blotted image")

    cr_mask = flag_cr(image.data, image.err, blot, snr1, snr2, scale1, scale2, backg, resample_data)

    # TODO is it necessary to do all this math for 1 possible log message?
    # sould a count of flagged pixels in the mask be a suitable replacement?
    # Count existing DO_NOT_USE pixels
    count_existing = np.count_nonzero(image.dq & DO_NOT_USE)

    # FIXME (or really "fixed") for a converted cube this used to overwrite
    # the dq "view" of the cube with a new dq array. This broke the link
    # between the converted ImageModel.dq and the original CubeModel.dq
    # which required extra work at the end of outlier detection to update
    # the cube. By modifying dq in-place we automatically update the cube.
    # I put this as a FIXME because this comment can be removed when
    # we don't convert cubes.
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    # Report number (and percent) of new DO_NOT_USE pixels found
    count_outlier = np.count_nonzero(image.dq & DO_NOT_USE)
    count_added = count_outlier - count_existing
    percent_cr = count_added / image.dq.size * 100
    log.info(f"New pixels flagged as outliers: {count_added} ({percent_cr:.2f}%)")


def _convert_inputs(inputs, good_bits, weight_type):
    """Convert input into datamodel required for processing.

    This base class works on imaging data, and relies on use of the
    ModelContainer class as the format needed for processing. However,
    the input may not always be a ModelContainer object, so this method
    will convert the input to a ModelContainer object for processing.
    Additionally, sub-classes may redefine this to set up the input as
    whatever format the sub-class needs for processing.

    """
    if isinstance(inputs, ModelContainer):
        return inputs
    input_models = ModelContainer()
    num_inputs = inputs.data.shape[0]
    log.debug("Converting CubeModel to ModelContainer with {} images".
              format(num_inputs))
    for i in range(inputs.data.shape[0]):
        image = datamodels.ImageModel(data=inputs.data[i],
                                      err=inputs.err[i],
                                      dq=inputs.dq[i])
        image.meta = inputs.meta
        image.wht = build_driz_weight(image,
                                      weight_type=weight_type,
                                      good_bits=good_bits)
        input_models.append(image)
    return input_models
