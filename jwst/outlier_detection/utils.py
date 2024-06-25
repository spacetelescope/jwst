"""
The ever-present utils sub-module. A home for all...
"""
import os
import warnings

import numpy as np
from astropy.stats import sigma_clip
from drizzle.cdrizzle import tblot
from scipy import ndimage

from jwst.resample.resample_utils import calc_gwcs_pixmap
from jwst.resample.resample import compute_image_pixel_area
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']


def _remove_file(fn):
    if isinstance(fn, str) and os.path.isfile(fn):
        os.remove(fn)
        log.info(f"Removing file {fn}")


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


def _absolute_subtract(array, tmp, out):
    tmp = np.abs(array - tmp)
    out = np.maximum(tmp, out)
    tmp = tmp * 0.
    return tmp, out


def flag_cr(sci_image, blot_data, snr="5.0 4.0", scale="1.2 0.7", backg=0,
            resample_data=True, **kwargs):
    """Masks outliers in science image by updating DQ in-place

    Mask blemishes in dithered data by comparing a science image
    with a model image and the derivative of the model image.

    Parameters
    ----------
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
    snr1, snr2 = [float(val) for val in snr.split()]
    scale1, scale2 = [float(val) for val in scale.split()]

    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (sci_image.meta.background.subtracted is False and
            sci_image.meta.background.level is not None):
        subtracted_background = sci_image.meta.background.level
        log.debug(f"Adding background level {subtracted_background} to blotted image")
    else:
        # No subtracted background.  Allow user-set value, which defaults to 0
        subtracted_background = backg

    sci_data = sci_image.data
    err_data = np.nan_to_num(sci_image.err)

    # create the outlier mask
    if resample_data:  # dithered outlier detection
        blot_deriv = _abs_deriv(blot_data)
        blot_data += subtracted_background
        diff_noise = np.abs(sci_data - blot_data)

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

    # Count existing DO_NOT_USE pixels
    count_existing = np.count_nonzero(sci_image.dq & DO_NOT_USE)

    # Update the DQ array in the input image.
    # FIXME (or really "fixed") for a converted cube this used to overwrite
    # the dq "view" of the cube with a new dq array. This broke the link
    # between the converted ImageModel.dq and the original CubeModel.dq
    # which required extra work at the end of outlier detection to update
    # the cube. By modifying dq in-place we automatically update the cube.
    # I put this as a FIXME because this comment can be removed when
    # we don't convert cubes.
    sci_image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    # Report number (and percent) of new DO_NOT_USE pixels found
    count_outlier = np.count_nonzero(sci_image.dq & DO_NOT_USE)
    count_added = count_outlier - count_existing
    percent_cr = count_added / sci_image.dq.size * 100
    log.info(f"New pixels flagged as outliers: {count_added} ({percent_cr:.2f}%)")


# FIXME (or fixed) interp and sinscl were "options" only when provided
# as part of the step spec (which becomes outlierpars). As neither was
# in the spec (and providing unknown arguments causes an error), these
# were never configurable and always defaulted to linear and 1.0
def gwcs_blot(median_model, blot_img):
    """
    Resample the output/resampled image to recreate an input image based on
    the input image's world coordinate system

    Parameters
    ----------
    median_model : `~stdatamodels.jwst.datamodels.JwstDataModel`

    blot_img : datamodel
        Datamodel containing header and WCS to define the 'blotted' image
    """
    blot_wcs = blot_img.meta.wcs

    # Compute the mapping between the input and output pixel coordinates
    pixmap = calc_gwcs_pixmap(blot_wcs, median_model.meta.wcs, blot_img.data.shape)
    log.debug("Pixmap shape: {}".format(pixmap[:, :, 0].shape))
    log.debug("Sci shape: {}".format(blot_img.data.shape))

    if 'SPECTRAL' not in blot_img.meta.wcs.output_frame.axes_type:
        input_pixflux_area = blot_img.meta.photometry.pixelarea_steradians
        input_pixel_area = compute_image_pixel_area(blot_img.meta.wcs)
        pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
    else:
        pix_ratio = 1.0
    log.info('Blotting {} <-- {}'.format(blot_img.data.shape, median_model.data.shape))

    outsci = np.zeros(blot_img.shape, dtype=np.float32)

    # Currently tblot cannot handle nans in the pixmap, so we need to give some
    # other value.  -1 is not optimal and may have side effects.  But this is
    # what we've been doing up until now, so more investigation is needed
    # before a change is made.  Preferably, fix tblot in drizzle.
    pixmap[np.isnan(pixmap)] = -1
    tblot(median_model.data, pixmap, outsci, scale=pix_ratio, kscale=1.0,
          interp='linear', exptime=1.0, misval=0.0, sinscl=1.0)

    return outsci


def detect_outliers(input_models, median_model, **kwargs):
    for image in input_models:
        if kwargs["resample_data"]:
            blot = gwcs_blot(median_model, image)
        else:
            blot = median_model.data
        # TODO save blot? are these actually useful?
        flag_cr(image, blot, **kwargs)
