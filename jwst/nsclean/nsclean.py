import logging
import numpy as np

import gwcs
from astropy.stats import sigma_clipped_stats

from jwst import datamodels
from jwst.assign_wcs import nirspec
from stdatamodels.jwst.datamodels import dqflags

from jwst.nsclean.lib import NSClean, NSCleanSubarray

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def mask_ifu_slices(input_model, mask):

    """Find pixels located within IFU slices and flag them in the
    mask, so that they do not get used.

    Parameters
    ----------
    input_model : data model object
        science data

    mask : 2D bool array
        input mask that will be updated

    Returns
    -------
    mask : 2D bool array
        output mask with additional flags for science pixels
    """

    log.info("Finding slice pixels for an IFU image")

    # Initialize global DQ map to all zero (OK to use)
    dqmap = np.zeros_like(input_model.dq)

    # Get the wcs objects for all IFU slices
    list_of_wcs = nirspec.nrs_ifu_wcs(input_model)

    # Loop over the IFU slices, finding the valid region for each
    for (k, ifu_wcs) in enumerate(list_of_wcs):

        # Construct array indexes for pixels in this slice
        x, y = gwcs.wcstools.grid_from_bounding_box(ifu_wcs.bounding_box,
                                                    step=(1, 1),
                                                    center=True)
        # Get the world coords for all pixels in this slice;
        # all we actually need are wavelengths
        coords = ifu_wcs(x, y)
        dq = dqmap[y.astype(int), x.astype(int)]
        wl = coords[2]
        # set non-NaN wavelength locations as do not use (one)
        valid = ~np.isnan(wl)
        dq = dq[valid]
        x = x[valid]
        y = y[valid]
        dq[:] = 1

        # Copy DQ for this slice into global DQ map
        dqmap[y.astype(int), x.astype(int)] = dq

    # Now set all non-zero locations in the mask to False (do not use)
    mask[dqmap==1] = False

    return mask


def mask_slits(input_model, mask):

    """Find pixels located within MOS or fixed slit footprints
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : data model object
        science data

    mask : 2D bool array
        input mask that will be updated

    Returns
    -------
    mask : 2D bool array
        output mask with additional flags for slit pixels
    """

    from jwst.extract_2d.nirspec import offset_wcs

    log.info("Finding slit/slitlet pixels")

    # Get the slit-to-msa frame transform from the WCS object
    slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')

    # Loop over the slits, marking all the pixels within each bounding
    # box as False (do not use) in the mask.
    # Note that for 3D masks (TSO mode), all planes will be set to the same value.
    for slit in slit2msa.slits:
        slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
        xlo, xhi, ylo, yhi = offset_wcs(slit_wcs)
        mask[..., ylo:yhi, xlo:xhi] = False

    return mask

    
def create_mask(input_model, mask_spectral_regions, n_sigma):
    """Create the pixel mask needed for setting which pixels to use
    for measuring 1/f noise.

    Parameters
    ----------
    input_model : data model object
        science data

    mask_spectral_regions : bool
        mask slit/slice regions defined in WCS

    n_sigma : float
        sigma threshold for masking outliers

    Returns
    -------
    mask : 2D or 3D bool array
        image mask

    nan_pix : array
        indexes of image locations with NaN values
    """
    exptype = input_model.meta.exposure.type.lower()

    # Initialize mask to all True. Subsequent operations will mask
    # out pixels that contain signal.
    # Note: mask will be 3D for BOTS mode data
    mask = np.full(np.shape(input_model.dq), True)

    # If IFU, mask all pixels contained in the IFU slices
    if exptype == 'nrs_ifu' and mask_spectral_regions:
        mask = mask_ifu_slices(input_model, mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if exptype in ['nrs_fixedslit', 'nrs_brightobj', 'nrs_msaspec'] and mask_spectral_regions:
        mask = mask_slits(input_model, mask)

    # If IFU or MOS, mask pixels affected by failed-open shutters
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        open_pix = input_model.dq & dqflags.pixel['MSA_FAILED_OPEN']
        mask[open_pix > 0] = False

    # Temporarily reset NaN pixels and mask them.
    # Save the list of NaN pixel coords, so that they can be reset at the end.
    nan_pix = np.isnan(input_model.data)
    input_model.data[nan_pix] = 0
    mask[nan_pix] = False

    # If IFU or MOS, mask the fixed-slit area of the image; uses hardwired indexes
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        mask[922:1116, :] = False

    # Use left/right reference pixel columns (first and last 4). Can only be
    # applied to data that uses all 2048 columns of the detector.
    if mask.shape[-1] == 2048:
        mask[..., :, :5] = True
        mask[..., :, -4:] = True

    # Mask outliers using sigma clipping stats.
    # For BOTS mode, which uses 3D data, loop over each integration separately.
    if len(input_model.data.shape) == 3:
        for i in range(input_model.data.shape[0]):
            _, median, sigma = sigma_clipped_stats(input_model.data[i], mask=~mask[i], mask_value=0, sigma=5.0)
            outliers = input_model.data[i] > (median + n_sigma * sigma)
            mask[i][outliers] = False
    else:
        _, median, sigma = sigma_clipped_stats(input_model.data, mask=~mask, mask_value=0, sigma=5.0)
        outliers = input_model.data > (median + n_sigma * sigma)
        mask[outliers] = False


    # Return the mask and the record of which pixels were NaN in the input;
    # it'll be needed later
    return mask, nan_pix


def clean_full_frame(detector, image, mask):
    """Clean a full-frame (2048x2048) image.

    Parameters
    ----------
    detector : str
        The name of the detector from which the data originate.

    image : 2D float array
        The image to be cleaned.

    mask : 2D bool array
        The mask that indicates which pixels are to be used in fitting.

    Returns
    -------
    cleaned_image : 2D float array
        The cleaned image.
    """

    # Instantiate the cleaner
    cleaner = NSClean(detector, mask)

    # Clean the image
    try:
        cleaned_image = cleaner.clean(image, buff=True)
    except np.linalg.LinAlgError:
        log.warning("Error cleaning image; step will be skipped")
        return None

    return cleaned_image


def clean_subarray(detector, image, mask):
    """Clean a subarray image.

    Parameters
    ----------
    detector : str
        The name of the detector from which the data originate.

    image : 2D float array
        The image to be cleaned.

    mask : 2D bool array
        The mask that indicates which pixels are to be used in fitting.

    Returns
    -------
    cleaned_image : 2D float array
        The cleaned image.
    """

    # Flip the image to detector coords. NRS1 requires a transpose
    # of the axes, while NRS2 requires a transpose and flip.
    if detector == "NRS1":
        image = image.transpose()
        mask = mask.transpose()
    else:
        image = image.transpose()[::-1]
        mask = mask.transpose()[::-1]

    # Instantiate the cleaner
    cleaner = NSCleanSubarray(image, mask)

    # Clean the image
    try:
        cleaned_image = cleaner.clean()
    except np.linalg.LinAlgError:
        log.warning("Error cleaning image; step will be skipped")
        return None

    # Restore the cleaned image to the science frame
    if detector == "NRS1":
        cleaned_image = cleaned_image.transpose()
    else:
        cleaned_image = cleaned_image[::-1].transpose()

    return cleaned_image


def do_correction(input_model, mask_spectral_regions, n_sigma, save_mask, user_mask):

    """Apply the NSClean 1/f noise correction

    Parameters
    ----------
    input_model : data model object
        science data to be corrected

    mask_spectral_regions : bool
        Mask slit/slice regions defined in WCS

    n_sigma : float
        n-sigma rejection level for finding outliers

    save_mask : bool
        switch to indicate whether the mask should be saved

    user_mask : str or None
        Path to user-supplied mask image

    Returns
    -------
    output_model : `~jwst.datamodel.JwstDataModel`
        corrected data

    mask_model : `~jwst.datamodel.JwstDataModel`
        pixel mask to be saved or None
    """

    detector = input_model.meta.instrument.detector.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')

    # Check for a valid input that we can work on
    if input_model.meta.subarray.name.upper() == "ALLSLITS":
        log.warning("Step cannot be applied to ALLSLITS subarray images")
        log.warning("Step will be skipped")
        input_model.meta.cal_step.nsclean = 'SKIPPED'
        return input_model, None

    output_model = input_model.copy()

    # Check for a user-supplied mask image. If so, use it.
    if user_mask is not None:
        mask_model = datamodels.open(user_mask)
        Mask = (mask_model.data.copy()).astype(np.bool_)

        # Reset and save list of NaN pixels in the input image
        nan_pix = np.isnan(input_model.data)
        input_model.data[nan_pix] = 0
        Mask[nan_pix] = False
    
    else:
        # Create the pixel mask that'll be used to indicate which pixels
        # to include in the 1/f noise measurements. Basically, we're setting
        # all illuminated pixels to False, so that they do not get used, and
        # setting all unilluminated pixels to True (so they DO get used).
        # For BOTS mode the mask will be 3D, to accommodate changes in masked
        # pixels per integration.
        log.info("Creating mask")
        Mask, nan_pix = create_mask(input_model, mask_spectral_regions, n_sigma)

        # Store the mask image in a model, if requested
        if save_mask:
            if len(Mask.shape) == 3:
                mask_model = datamodels.CubeModel(data=Mask)
            else:
                mask_model = datamodels.ImageModel(data=Mask)
        else:
            mask_model = None

    log.info(f"Cleaning image {input_model.meta.filename}")

    # Setup for handling 2D or 3D inputs
    if len(input_model.data.shape) == 3:
        nints = input_model.data.shape[0]
        # Check for 3D mask
        if len(Mask.shape) == 2:
            log.warning("Data are 3D, but mask is 2D. Step will be skipped.")
            output_model.meta.cal_step.nsclean = 'SKIPPED'
            return output_model, None
    else:
        nints = 1

    # Loop over integrations (even if there's only 1)
    for i in range(nints):
        log.debug(f" working on integration {i+1}")
        if len(input_model.data.shape) == 3:
            image = np.float32(input_model.data[i])
            mask = Mask[i]
        else:
            image = np.float32(input_model.data)
            mask = Mask

        if input_model.data.shape[-2:] == (2048, 2048):
            # Clean a full-frame image
            cleaned_image = clean_full_frame(detector, image, mask)
        else:
            # Clean a subarray image
            cleaned_image = clean_subarray(detector, image, mask)

        # Check for failure
        if cleaned_image is None:
            output_model.meta.cal_step.nsclean = 'SKIPPED'
            break
        else:
            # Store the cleaned image in the output model
            if len(output_model.data.shape) == 3:
                output_model.data[i] = cleaned_image
            else:
                output_model.data = cleaned_image

    # Restore NaN's from original image
    output_model.data[nan_pix] = np.nan

    # Set completion status
    output_model.meta.cal_step.nsclean = 'COMPLETE'

    return output_model, mask_model
