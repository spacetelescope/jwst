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

    sigma: float, optional
        Number of standard deviations to use for both the lower
        and upper clipping limits.

    maxiters: int or None, optional
        Maximum number of sigma-clipping iterations to perform

    bkg_model: JWST data model
        background data model

    Returns
    -------
    result: JWST data model
        background-subtracted target data model

    """
    # Subtract the average background from the member
    log.debug(' subtracting avg bkg from {}'.format(input_model.meta.filename))

    # Compute the average of the background images associated with
    # the target exposure
    if input_model.meta.model_type == "ImageModel":
        data_dim = 2  # rate
    elif input_model.meta.model_type == "CubeModel":
        data_dim = 3  # rateints
    else:
        log.error("Input target exposure must be either a CubeModel or an ImageModel.")
        raise ValueError("Expected either a CubeModel or an ImageModel.")

    bkg_model = average_background(bkg_list, sigma, maxiters, data_dim)

    if data_dim == 2:  # rate
        result = subtract_images.subtract(input_model, bkg_model)
    elif data_dim == 3:  # rateints
        result = subtract_images.subtract_ints(input_model, bkg_model)

    # We're done. Return the background model and result.
    return bkg_model, result


def average_background(bkg_list, sigma, maxiters, data_dim):
    """
    Short Summary
    -------------
    Average multiple background exposures into a combined data model.

    Long Summary
    ------------
    For 2D (rate) background exposures, the average background will be calculated
    as the mean of the sigma-clipped data in the background expoosures.

    For 3D (rateint) background exposures, each exposure can have an arbitrary
    number of integrations. First the average (over integrations) of the clipped
    mean will be calculated for each exposure. Then the mean of these average
    backgrounds will be calculated as the final average background image. This
    final average will be subtracted from each integration of the SCI exposure.

    Parameters:
    -----------
    bkg_list: filename list
        List of background exposure file names

    sigma: float, optional
        Number of standard deviations to use for both the lower
        and upper clipping limits.

    maxiters: int or None, optional
        Maximum number of sigma-clipping iterations to perform

    data_dim: int
        Dimension of data arrays in input files. For input rate files: 2.
        For input rateint files: 3,

    Returns:
    --------
    avg_bkg: data model
        The averaged background exposure
    """
    num_bkg = len(bkg_list)
    avg_bkg = None

    if data_dim == 2:  # background rate input
        cdata = None
        # Loop over the images to be used as background
        for i, bkg_file in enumerate(bkg_list):
            log.info(f'Accumulate bkg from {bkg_file}')
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
        log.debug('clip with sigma={} maxiters={}'.format(sigma, maxiters))
        mdata = sigma_clip(cdata, sigma=sigma, maxiters=maxiters, axis=0)

        # Compute the mean of the non-clipped values
        avg_bkg.data = mdata.mean(axis=0).data

        # Mask the ERR values using the data mask
        merr = np.ma.masked_array(cerr, mask=mdata.mask)

        # Compute the combined ERR as the uncertainty in the mean
        avg_bkg.err = (np.sqrt(merr.sum(axis=0)) / (num_bkg - merr.mask.sum(axis=0))).data

        return avg_bkg

    else:  # data_dim = 3; background rateints input
        temp_avg_bkg = None  # accumulating Cube model; file is along axis 0

        # Loop over the images to be used as background
        for i, bkg_file in enumerate(bkg_list):
            log.info(f'Accumulate bkg from {bkg_file}')
            bkg_model = datamodels.CubeModel(bkg_file)  # axis=0 is for integrations

            # Initialize a bunch of stuff if necessary
            if temp_avg_bkg is None:
                image_shape = bkg_model.data.shape[1:]
                avg_bkg = datamodels.ImageModel(image_shape)
                accum_dq_arr = np.zeros((image_shape), dtype=np.int32)
                temp_avg_bkg = datamodels.CubeModel((num_bkg,) + image_shape)

            # Sigma clip the bkg model's data along the integration axis
            sc_bkg_data = sigma_clip(bkg_model.data, sigma=sigma, maxiters=maxiters, axis=0)

            # Accumulate the integ-averaged clipped data in the temp Cubemodel for file
            temp_avg_bkg.data[i, :, :] = sc_bkg_data.mean(axis=0)

            # Collapse the DQ by doing a bitwise_OR over all integrations in the file
            for i_nint in range(bkg_model.dq.shape[0]):
                accum_dq_arr = np.bitwise_or(bkg_model.dq[i_nint, :, :], accum_dq_arr)

            # Accumulate the collapsed DQ in the temp Cubemodel for file
            temp_avg_bkg.dq[i, :, :] = accum_dq_arr

        # Average data averages over all files
        avg_bkg.data = temp_avg_bkg.data.mean(axis=0)  # mean along file axis

        # Do bitwise_OR over all of the files to give 2D dq array
        temp_dq_arr = temp_avg_bkg.dq[0, :, :]
        for i in range(num_bkg):
            temp_dq_arr = np.bitwise_or(temp_avg_bkg.dq[i, :, :], temp_dq_arr)

        avg_bkg.dq = temp_dq_arr

        return avg_bkg


def subtract_wfss_bkg(input_model, bkg_filename, wl_range_name, mmag_extract=None):
    """
    Short Summary
    -------------
    Scale and subtract a background reference image from WFSS/GRISM data.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    bkg_filename: str
        name of master background file for WFSS/GRISM

    wl_range_name: str
        name of wavelengthrange reference file

    mmag_extract: float
        minimum abmag of grism objects to extract

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
        bkg_mask = mask_from_source_cat(input_model, wl_range_name, mmag_extract)
        if bkg_mask.sum() < 100:
            log.warning("Not enough background pixels to work with.")
            log.warning("Step will be SKIPPED.")
            return None
    else:
        bkg_mask = np.ones(input_model.data.shape, dtype=bool)

    # Compute the mean values of science image and background reference
    # image, including only regions where there are no identified sources.
    # Exclude pixel values in the lower and upper 25% of the histogram.
    lowlim = 25.
    highlim = 75.
    sci_mean = robust_mean(input_model.data[bkg_mask],
                           lowlim=lowlim, highlim=highlim)
    bkg_mean = robust_mean(bkg_ref.data[bkg_mask],
                           lowlim=lowlim, highlim=highlim)

    log.debug("mean of [{}, {}] percentile grism image = {}"
              .format(lowlim, highlim, sci_mean))
    log.debug("mean of [{}, {}] percentile background image = {}"
              .format(lowlim, highlim, bkg_mean))

    result = input_model.copy()
    if bkg_mean != 0.:
        subtract_this = (sci_mean / bkg_mean) * bkg_ref.data
        result.data = input_model.data - subtract_this
        log.info(f"Average of background image subtracted = {subtract_this.mean(dtype=float)}")
    else:
        log.warning("Background image has zero mean; nothing will be subtracted.")
    result.dq = np.bitwise_or(input_model.dq, bkg_ref.dq)

    bkg_ref.close()

    return result


def no_NaN(model, fill_value=0.):
    """
    Short Summary
    -------------
    Replace NaNs with a harmless value.

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


def mask_from_source_cat(input_model, wl_range_name, mmag_extract=None):
    """
    Short Summary
    -------------
    Create a mask that is False within bounding boxes of sources.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    wl_range_name: str
        Name of the wavelengthrange reference file

    mmag_extract: float
        Minimum abmag of grism objects to extract

    Returns
    -------
    bkg_mask: ndarray
        Boolean mask:  True for background, False for pixels that are
        inside at least one of the source regions defined in the source
        catalog.
    """

    shape = input_model.data.shape
    bkg_mask = np.ones(shape, dtype=bool)

    reference_files = {"wavelengthrange": wl_range_name}
    grism_obj_list = create_grism_bbox(input_model, reference_files, mmag_extract)

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
    """
    Short Summary
    -------------
    Compute a mean value, excluding outliers.

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
    mean_value = x[mask].mean(dtype=float)

    return mean_value
