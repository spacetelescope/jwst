import math
import numpy as np

from .. import datamodels
from . import subtract_images
from ..assign_wcs.util import create_grism_bbox
from astropy.stats import sigma_clip

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def background_sub(input_model, bkg_list, sigma, maxiters):

    """
    Short Summary
    -------------
    Subtract the background signal from a JWST exposure by subtracting
    the average of one or more background exposures from it.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    bkg_list: filename list
        list of background exposure file names

    Returns
    -------
    result: JWST data model
        background-subtracted target data model

    """

    # Compute the average of the background images associated with
    # the target exposure
    bkg_model = average_background(bkg_list, sigma, maxiters)

    # Subtract the average background from the member
    log.debug(' subtracting avg bkg from {}'.format(input_model.meta.filename))
    result = subtract_images.subtract(input_model, bkg_model)

    # Close the average background image and update the step status
    bkg_model.close()

    # We're done. Return the result.
    return result


def average_background(bkg_list, sigma, maxiters):

    """
    Average multiple background exposures into a combined data model

    Parameters:
    -----------

    bkg_list: filename list
        List of background exposure file names

    Returns:
    --------

    avg_bkg: data model
        The averaged background exposure

    """

    num_bkg = len(bkg_list)
    avg_bkg = None
    cdata = None

    # Loop over the images to be used as background
    for i, bkg_file in enumerate(bkg_list):
        log.debug(' Accumulate bkg from {}'.format(bkg_file))
        bkg_model = datamodels.ImageModel(bkg_file)

        # Initialize the avg_bkg model, if necessary
        if avg_bkg is None:
            avg_bkg = datamodels.ImageModel(bkg_model.shape)

        if cdata is None:
            cdata = np.zeros(((num_bkg,) + bkg_model.shape))
            cerr = cdata.copy()

        # Accumulate the data from this background image
        cdata[i] = bkg_model.data
        cerr[i] = bkg_model.err * bkg_model.err
        avg_bkg.dq = np.bitwise_or(avg_bkg.dq, bkg_model.dq)

        bkg_model.close()

    # Clip the background data
    log.debug(' clip with sigma={} maxiters={}'.format(sigma, maxiters))
    mdata = sigma_clip(cdata, sigma=sigma, maxiters=maxiters, axis=0)

    # Compute the mean of the non-clipped values
    avg_bkg.data = mdata.mean(axis=0).data

    # Mask the ERR values using the data mask
    merr = np.ma.masked_array(cerr, mask=mdata.mask)

    # Compute the combined ERR as the uncertainty in the mean
    avg_bkg.err = np.sqrt(merr.sum(axis=0)) / (num_bkg - merr.mask.sum(axis=0))

    return avg_bkg


def subtract_wfss_bkg(input_model, bkg_filename, wl_range_name):
    """Scale and subtract a background reference image from WFSS/GRISM data.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    bkg_filename: str
        name of master background file for WFSS/GRISM

    wl_range_name: str
        name of wavelengthrange reference file

    Returns
    -------
    result: JWST data model
        background-subtracted target data model
    """

    bkg_ref = datamodels.open(bkg_filename)

    if hasattr(input_model.meta, "source_catalog"):
        got_catalog = True
    else:
        log.warning("No source_catalog found in input.meta.")
        got_catalog = False

    # If there are NaNs, we have to replace them with something harmless.
    bkg_ref = no_NaN(bkg_ref)

    # Create a mask from the source catalog, True where there are no sources,
    # i.e. in regions we can use as background.
    if got_catalog:
        bkg_mask = mask_from_source_cat(input_model, wl_range_name)
    else:
        bkg_mask = np.ones(input_model.data.shape, dtype=np.bool)

    # Compute the mean values of science image and background reference
    # image, including only regions where there are no identified sources.
    # Exclude pixel values in the lower and upper 25% of the histogram.
    lowlim = 25.
    highlim = 75.
    sci_mean = robust_mean(input_model.data[bkg_mask],
                           lowlim=lowlim, highlim=highlim)
    bkg_mean = robust_mean(bkg_ref.data[bkg_mask],
                           lowlim=lowlim, highlim=highlim)

    log.debug("mean of [{}, {}] percentile science data = {}"
              .format(lowlim, highlim, sci_mean))
    log.debug("mean of [{}, {}] percentile background image = {}"
              .format(lowlim, highlim, bkg_mean))

    result = input_model.copy()
    if bkg_mean != 0.:
        subtract_this = (sci_mean / bkg_mean) * bkg_ref.data
        result.data = input_model.data - subtract_this
        log.debug("Average of values subtracted = {}"
                  .format(subtract_this.mean(dtype=np.float)))
    else:
        log.warning("Background file has zero mean; "
                    "nothing will be subtracted.")
    result.dq = np.bitwise_or(input_model.dq, bkg_ref.dq)

    bkg_ref.close()

    return result


def no_NaN(model, fill_value=0.):
    """Replace NaNs with a harmless value.

    Parameters
    ----------
    model: JWST data model
        Reference file model.

    fill_value: float
        NaNs will be replaced with this value.

    Returns
    -------
    result: JWST data model
        Reference file model without NaNs in data array.
    """

    mask = np.isnan(model.data)
    if mask.sum(dtype=np.intp) == 0:
        return model
    else:
        temp = model.copy()
        temp.data[mask] = fill_value
        return temp


def mask_from_source_cat(input_model, wl_range_name):
    """Create a mask that is False within bounding boxes of sources.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    wl_range_name: str
        Name of the wavelengthrange reference file

    Returns
    -------
    bkg_mask: ndarray
        Boolean mask:  True for background, False for pixels that are
        inside at least one of the source regions defined in the source
        catalog.
    """

    shape = input_model.data.shape
    bkg_mask = np.ones(shape, dtype=np.bool)

    reference_files = {"wavelengthrange": wl_range_name}
    grism_obj_list = create_grism_bbox(input_model, reference_files)

    for obj in grism_obj_list:
        order_bounding = obj.order_bounding
        for order in order_bounding.keys():
            ((ymin, ymax), (xmin, xmax)) = order_bounding[order]
            xmin = int(math.floor(xmin))
            xmax = int(math.ceil(xmax)) + 1     # convert to slice limit
            ymin = int(math.floor(ymin))
            ymax = int(math.ceil(ymax)) + 1
            xmin = max(xmin, 0)
            xmax = min(xmax, shape[-1])
            ymin = max(ymin, 0)
            ymax = min(ymax, shape[-2])
            bkg_mask[..., ymin:ymax, xmin:xmax] = False

    return bkg_mask


def robust_mean(x, lowlim=25., highlim=75.):
    """Compute a mean value, excluding outliers.

    Parameters
    ----------
    x: ndarray
        The array for which we want a mean value.

    lowlim: float
        The lower `lowlim` percent of the data will not be used when
        computing the mean.

    highlim: float
        The upper `highlim` percent of the data will not be used when
        computing the mean.

    Returns
    -------
    mean_value: float
        The mean of `x`, excluding data outside `lowlim` to `highlim`
        percentile limits.
    """

    limits = np.percentile(x, (lowlim, highlim))
    mask = np.logical_and(x >= limits[0], x <= limits[1])
    mean_value = x[mask].mean(dtype=np.float)

    return mean_value
