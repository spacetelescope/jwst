import copy
import numpy as np
import warnings

from stdatamodels.jwst import datamodels

from . import subtract_images
from astropy.stats import sigma_clip
from astropy.utils.exceptions import AstropyUserWarning

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ImageSubsetArray:
    """
    Keep track of where different images relate to each other.

    This includes using e.g. subarray observations.
    """

    def __init__(self, model):
        """
        Initialize the class.

        Parameters
        ----------
        model : ImageModel
            Input datamodel
        """
        im = datamodels.open(model)

        # Make sure xstart/ystart/xsize/ysize exist in the model metadata
        if (
            im.meta.subarray.xstart is None
            or im.meta.subarray.ystart is None
            or im.meta.subarray.xsize is None
            or im.meta.subarray.ysize is None
        ):
            # If the ref file is full-frame, set the missing params to
            # default values
            if im.meta.instrument.name.upper() == "MIRI":
                if im.data.shape[-1] == 1032 and im.data.shape[-2] == 1024:
                    im.meta.subarray.xstart = 1
                    im.meta.subarray.ystart = 1
                    im.meta.subarray.xsize = 1032
                    im.meta.subarray.ysize = 1024
                else:
                    raise ValueError("xstart or ystart metadata values not found in model")
            else:
                if im.data.shape[-1] == 2048 and im.data.shape[-2] == 2048:
                    im.meta.subarray.xstart = 1
                    im.meta.subarray.ystart = 1
                    im.meta.subarray.xsize = 2048
                    im.meta.subarray.ysize = 2048
                else:
                    raise ValueError("xstart or ystart metadata values not found in model")

        # Account for the fact these things are 1-indexed
        self.jmin = im.meta.subarray.ystart - 1
        self.jmax = im.meta.subarray.ysize + self.jmin
        self.imin = im.meta.subarray.xstart - 1
        self.imax = im.meta.subarray.xsize + self.imin

        # Check the dimensionality of this thing
        self.im_dim = len(im.data.shape)

        self.data = im.data
        self.err = im.err
        self.dq = im.dq

        im.close()

    def overlaps(self, other):
        """
        Find whether this subset and another overlap.

        Parameters
        ----------
        other : ImageModel
            Input model

        Returns
        -------
        bool
            Return the boolean opposite of the condition met
        """
        return not (
            self.imax <= other.imin
            or other.imax <= self.imin
            or self.jmax <= other.jmin
            or other.jmax <= self.jmin
        )

    def get_subset_array(self, other):
        """
        Pull out the overlapping SCI, ERR, and DQ data for two models.

        Parameters
        ----------
        other : ImageModel
            Input model

        Returns
        -------
        data_overlap : ndarray
            Data overlap of the science arrays

        err_overlap : ndarray
            Overlap of the error arrays

        dq_overlap : ndarray
            Overlap of the DQ arrays
        """
        imin = max(self.imin, other.imin)
        imax = min(self.imax, other.imax)
        jmin = max(self.jmin, other.jmin)
        jmax = min(self.jmax, other.jmax)

        if imax < imin:
            imax = imin

        if jmax < jmin:
            jmax = jmin

        # To ensure that we can mix and match subarray obs, take
        # the x/y shape from self.data. To ensure we can work
        # with mismatched NINTS, if we have that third dimension
        # use the value from other
        data_shape = list(self.data.shape)
        if self.im_dim == 3:
            data_shape[0] = other.data.shape[0]

        # Set up arrays, NaN out data/err for sigma clipping, keep DQ as 0 for bitwise_or
        data_overlap = np.full(data_shape, np.nan, dtype=other.data.dtype)
        err_overlap = np.full(data_shape, np.nan, dtype=other.data.dtype)
        dq_overlap = np.zeros(data_shape, dtype=np.uint32)

        if self.im_dim == 2:
            idx = (
                slice(jmin - other.jmin, jmax - other.jmin),
                slice(imin - other.imin, imax - other.imin),
            )
        else:
            idx = (
                slice(None, None),
                slice(jmin - other.jmin, jmax - other.jmin),
                slice(imin - other.imin, imax - other.imin),
            )

        # Pull the values from the background obs
        data_cutout = other.data[idx]
        err_cutout = other.err[idx]
        dq_cutout = other.dq[idx]

        # Put them into the right place in the original image array shape
        data_overlap[..., : data_cutout.shape[-2], : data_cutout.shape[-1]] = copy.deepcopy(
            data_cutout
        )
        err_overlap[..., : data_cutout.shape[-2], : data_cutout.shape[-1]] = copy.deepcopy(
            err_cutout
        )
        dq_overlap[..., : data_cutout.shape[-2], : data_cutout.shape[-1]] = copy.deepcopy(dq_cutout)

        return data_overlap, err_overlap, dq_overlap


def background_sub(input_model, bkg_list, sigma, maxiters):
    """
    Subtract background signal from an exposure.

    This will be done by subtracting the average of one or more
    background exposures from it.

    Parameters
    ----------
    input_model : ImageModel
        Input target exposure data model

    bkg_list : list
        Filename list of background exposures

    sigma : float
        Number of standard deviations to use for both the lower
        and upper clipping limits.

    maxiters : int or None
        Maximum number of sigma-clipping iterations to perform

    Returns
    -------
    bkg_model : ImageModel
        Background data model

    result : ImageModel
        Background-subtracted target data model
    """
    # Compute the average of the background images associated with
    # the target exposure
    bkg_model = average_background(
        input_model,
        bkg_list,
        sigma,
        maxiters,
    )

    # Subtract the average background from the member
    log.info(f"Subtracting avg bkg from {input_model.meta.filename}")

    result = subtract_images.subtract(input_model, bkg_model)

    # We're done. Return the background model and the result.
    return bkg_model, result


def average_background(input_model, bkg_list, sigma, maxiters):
    """
    Average multiple background exposures into a combined data model.

    Processes backgrounds from various DataModel types, including those
    having 2D (rate) or 3D (rateints) backgrounds.

    Parameters
    ----------
    input_model : ImageModel
        Input target exposure data model

    bkg_list : list
        File name list of background exposures

    sigma : float, optional
        Number of standard deviations to use for both the lower
        and upper clipping limits.

    maxiters : int or None, optional
        Maximum number of sigma-clipping iterations to perform

    Returns
    -------
    avg_bkg : ImageModel
        The averaged background exposure
    """
    # Determine the dimensionality of the input file
    bkg_dim = len(input_model.data.shape)
    image_shape = input_model.data.shape[-2:]

    im_array = ImageSubsetArray(input_model)

    avg_bkg = datamodels.ImageModel(image_shape)
    num_bkg = len(bkg_list)
    cdata = np.zeros((num_bkg,) + image_shape)
    cerr = cdata.copy()

    if bkg_dim == 3:
        accum_dq_arr = np.zeros((image_shape), dtype=np.uint32)

    # Loop over the images to be used as background
    for i, bkg_file in enumerate(bkg_list):
        log.info(f"Accumulate bkg from {bkg_file}")

        bkg_array = ImageSubsetArray(bkg_file)

        if not bkg_array.overlaps(im_array):
            # We don't overlap, so put in a bunch of NaNs so sigma-clip
            # isn't affected and move on
            log.debug(f"{bkg_file} does not overlap input image")
            cdata[i] = np.ones(image_shape) * np.nan
            cerr[i] = np.ones(image_shape) * np.nan
            continue

        bkg_data, bkg_err, bkg_dq = im_array.get_subset_array(bkg_array)

        if bkg_dim == 2:
            # Accumulate the data from this background image
            cdata[i] = bkg_data  # 2D slice
            cerr[i] = bkg_err * bkg_err
            avg_bkg.dq = np.bitwise_or(avg_bkg.dq, bkg_dq)

        if bkg_dim == 3:
            # Sigma clip the bkg model's data and err along the integration axis
            with warnings.catch_warnings():
                # clipping NaNs and infs is the expected behavior
                warnings.filterwarnings(
                    "ignore", category=AstropyUserWarning, message=".*automatically clipped.*"
                )
                sc_bkg_data = sigma_clip(bkg_data, sigma=sigma, maxiters=maxiters, axis=0)
                sc_bkg_err = sigma_clip(bkg_err * bkg_err, sigma=sigma, maxiters=maxiters, axis=0)

            # Accumulate the integ-averaged clipped data and err for the file
            cdata[i] = sc_bkg_data.mean(axis=0)
            cerr[i] = sc_bkg_err.mean(axis=0)

            # Collapse the DQ by doing a bitwise_OR over all integrations
            for i_nint in range(bkg_dq.shape[0]):
                accum_dq_arr = np.bitwise_or(bkg_dq[i_nint, :, :], accum_dq_arr)
            avg_bkg.dq = np.bitwise_or(avg_bkg.dq, accum_dq_arr)

    # Clip the background data
    log.debug(f"clip with sigma={sigma} maxiters={maxiters}")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=AstropyUserWarning, message=".*automatically clipped.*"
        )
        mdata = sigma_clip(cdata, sigma=sigma, maxiters=maxiters, axis=0)

    # Compute the mean of the non-clipped values
    avg_bkg.data = mdata.mean(axis=0).data

    # Mask the ERR values using the data mask
    merr = np.ma.masked_array(cerr, mask=mdata.mask)

    # Compute the combined ERR as the uncertainty in the mean
    avg_bkg.err = (np.sqrt(merr.sum(axis=0)) / (num_bkg - merr.mask.sum(axis=0))).data

    return avg_bkg
