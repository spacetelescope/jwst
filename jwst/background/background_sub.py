from __future__ import division

import math
import numpy as np
from astropy.table import Table

from .. import datamodels
from . import subtract_images

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def background_sub(input_model, bkg_list):

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

    if input_model.meta.exposure.type in ["NIS_WFSS", "NRC_GRISM"]:

        # This is a reference file.
        if len(bkg_list) > 0:
            bkg_ref = datamodels.open(bkg_list[0])
            result = subtract_wfss_bkg(input_model, bkg_ref)
            bkg_ref.close()
        else:
            log.info("No reference file available; "
                     "skipping background subtraction.")
            result = input_model.copy()
            result.meta.cal_step.back_sub = 'SKIPPED'
            return result

    else:

        # Compute the average of the background images associated with
        # the target exposure
        bkg_model = average_background(bkg_list)

        # Subtract the average background from the member
        log.debug(' subtracting avg bkg from %s', input_model.meta.filename)
        result = subtract_images.subtract(input_model, bkg_model)

        # Close the average background image and update the step status
        bkg_model.close()

    # We're done. Return the result.
    return result


def average_background(bkg_list):

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

    avg_bkg = None

    # Loop over the images to be used as background
    for bkg_file in bkg_list:
        log.debug(' Accumulate bkg from %s', bkg_file)
        bkg_model = datamodels.ImageModel(bkg_file)

        # Initialize the avg_bkg model, if necessary
        if avg_bkg is None:
            avg_bkg = datamodels.ImageModel(bkg_model.data.shape)

        # Accumulate the data from this background image
        avg_bkg.data += bkg_model.data
        avg_bkg.err += bkg_model.err * bkg_model.err
        avg_bkg.dq = np.bitwise_or(avg_bkg.dq, bkg_model.dq)

        bkg_model.close()

    # Average the data in the accumulated background image
    num_bkg = len(bkg_list)
    avg_bkg.data = avg_bkg.data / num_bkg  # sci is normal average
    avg_bkg.err = np.sqrt(avg_bkg.err) / num_bkg  # err is uncertainty in the mean

    return avg_bkg


def subtract_wfss_bkg(input_model, bkg_ref):
    """Scale and subtract a background reference image from WFSS/GRISM data.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    Returns
    -------
    result: JWST data model
        background-subtracted target data model
    """

    if hasattr(input_model.meta, "source_catalog"):
        source_catalog = input_model.meta.source_catalog.filename
        got_catalog = True
    else:
        log.warning("No source_catalog found in input.meta.")
        got_catalog = False

    # If there are NaNs, we have to replace them with something harmless.
    bkg_ref = no_NaN(bkg_ref)

    # Create a mask from the source catalog, True where there are no sources,
    # i.e. in regions we can use as background.
    if got_catalog:
        bkg_mask = mask_from_source_cat(input_model.data.shape, source_catalog)
    else:
        bkg_mask = np.ones(input_model.data.shape, dtype=np.bool)

    # Compute the mean values of science image and background reference
    # image, including only regions where there are no identified sources.
    # Exclude pixel values in the lower 25% and upper 25% of the histogram.
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
        result.data = input_model.data - (sci_mean / bkg_mean) * bkg_ref.data
    else:
        log.warning("Background file has zero mean; "
                    "nothing will be subtracted.")
    result.dq = np.bitwise_or(input_model.dq, bkg_ref.dq)
    return result


def no_NaN(model, fill_value=0.):
    """Replace NaN's with a harmless value.

    Parameters
    ----------
    model: JWST data model
        Reference file model.

    Returns
    -------
    result: JWST data model
        Reference file model without NaN's in data array.
    """

    mask = (np.isnan(model.data))
    temp = model.copy()
    temp.data[mask] = fill_value
    return temp


def mask_from_source_cat(shape, source_catalog):
    """Create a mask that is False within bounding boxes of sources.

    Parameters
    ----------
    shape: tuple
        The boolean mask will be created with this shape.

    source_catalog: str
        Name of the source catalog (a .ecsv file).

    Returns
    -------
    bkg_mask: ndarray
        Boolean mask:  True for background, False for pixels that are
        inside at least one of the source regions defined in the source
        catalog.
    """

    # This was mostly copied from get_object_info() in assign_wcs/util.py

    bkg_mask = np.ones(shape, dtype=np.bool)

    catalog = Table.read(source_catalog, format="ascii.ecsv")
    if len(catalog) < 1:
        log.error("The source catalog is empty.")
        return bkg_mask
    else:
        msg = "Missing keys in catalog file: "
        ok = True
        row = catalog[0]
        try:
            xmin = row['xmin']
        except KeyError:
            msg = msg + " xmin"
            ok = False
        try:
            xmax = row['xmax']
        except KeyError:
            msg = msg + " xmax"
            ok = False
        try:
            ymin = row['ymin']
        except KeyError:
            msg = msg + " ymin"
            ok = False
        try:
            ymax = row['ymax']
        except KeyError:
            msg = msg + " ymax"
            ok = False
        if not ok:
            log.error(msg)
            return bkg_mask

    for row in catalog:
        # Read values, and convert from inclusive limits to slice limits.
        xmin = int(math.floor(row['xmin']))
        xmax = int(math.ceil(row['xmax'])) + 1
        ymin = int(math.floor(row['ymin']))
        ymax = int(math.ceil(row['ymax'])) + 1
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
