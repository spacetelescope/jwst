import logging
import numpy as np

from astropy.stats import sigma_clipped_stats

from jwst import datamodels
from .. assign_wcs import nirspec
from stdatamodels.jwst.datamodels import dqflags

from jwst.nsclean.lib import NSClean

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def mask_ifu_slices(input_model, Mask):

    """Find pixels located within IFU slices and flag them in the
    mask, so that they do not get used.

    Parameters
    ----------
    input_model : data model object
        science data

    Mask : 2D boolean array
        input mask that will be updated

    Returns
    -------
    Mask : 2D boolean array
        output mask with additional flags for science pixels
    """

    import gwcs

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
    Mask[dqmap==1] = False

    return Mask


def mask_slits(input_model, Mask):

    """Find pixels located within MOS or fixed slit footprints
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : data model object
        science data

    Mask : 2D boolean array
        input mask that will be updated

    Returns
    -------
    Mask : 2D boolean array
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
        Mask[..., ylo:yhi, xlo:xhi] = False

    return Mask

    
def create_mask(input_model, n_sigma):
    """Create the pixel mask needed for setting which pixels to use
    for measuring 1/f noise.

    Parameters
    ----------
    input_model : data model object
        science data

    n_sigma : float
        sigma threshold for masking outliers

    Returns
    -------
    Mask : 2D boolean array
        image mask

    nan_pix : array
        indexes of image locations with NaN values
    """
    exptype = input_model.meta.exposure.type.lower()

    # Initialize mask to all True. Subsequent operations will mask
    # out pixels that contain signal.
    # Note: Mask will be 3D for BOTS mode data
    Mask = np.full(np.shape(input_model.dq), True)

    # If IFU, mask all pixels contained in the IFU slices
    if exptype == 'nrs_ifu':
        Mask = mask_ifu_slices(input_model, Mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if exptype in ['nrs_fixedslit', 'nrs_brightobj', 'nrs_msaspec']:
        Mask = mask_slits(input_model, Mask)

    # If IFU or MOS, mask pixels affected by failed-open shutters
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        open_pix = np.where(input_model.dq & dqflags.pixel['MSA_FAILED_OPEN'])
        Mask[open_pix] = False

    # Temporarily reset NaN pixels and mask them.
    # Save the list of NaN pixel coords, so that they can be reset at the end.
    nan_pix = np.where(np.isnan(input_model.data))
    input_model.data[nan_pix] = 0
    Mask[nan_pix] = False

    # If IFU or MOS, mask the fixed-slit area of the image; uses hardwired indexes
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        Mask[922:1116, :] = False

    # Use left/right reference pixel columns (first and last 4). Can only be
    # applied to data that uses all 2048 columns of the detector.
    if Mask.shape[-1] == 2048:
        Mask[..., :, :5] = True
        Mask[..., :, -4:] = True

    # Mask outliers using sigma clipping stats.
    # For BOTS mode, which uses 3D data, loop over each integration separately.
    if len(input_model.data.shape) == 3:
        for i in range(input_model.data.shape[0]):
            _, median, sigma = sigma_clipped_stats(input_model.data[i], mask=~Mask[i], mask_value=0, sigma=5.0)
            outliers = np.where(input_model.data[i] > (median + n_sigma * sigma))
            Mask[i][outliers] = False
    else:
        _, median, sigma = sigma_clipped_stats(input_model.data, mask=~Mask, mask_value=0, sigma=5.0)
        outliers = np.where(input_model.data > (median + n_sigma * sigma))
        Mask[outliers] = False


    # Return the mask and the record of which pixels were NaN in the input;
    # it'll be needed later
    return Mask, nan_pix


def do_correction(input_model, n_sigma, save_mask, user_mask):

    """Apply the NSClean 1/f noise correction

    Parameters
    ----------
    input_model : data model object
        science data to be corrected

    n_sigma : float
        n-sigma rejection level for finding outliers

    save_mask : boolean
        switch to indicate whether the mask should be saved

    user_mask : string or None
        Path to user-supplied mask image

    Returns
    -------
    output_model : jwst.datamodel.JwstDataModel
        corrected data

    mask_model : jwst.datamodel.JwstDataModel
        pixel mask to be saved or None
    """

    detector = input_model.meta.instrument.detector.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')
    output_model = input_model.copy()

    # Check for full readout length in dispersion direction
    if input_model.data.shape[-1] != 2048:
        log.warning("Can't clean subarrays of this size; step will be skipped")
        output_model.meta.cal_step.nsclean = 'SKIPPED'
        return output_model

    # Check for a user-supplied mask image. If so, use it.
    if user_mask is not None:
        mask_model = datamodels.open(user_mask)
        Mask = mask_model.data.copy()
    
    else:
        # Create the pixel mask that'll be used to indicate which pixels
        # to include in the 1/f noise measurements. Basically, we're setting
        # all illuminated pixels to False, so that they do not get used, and
        # setting all unilluminated pixels to True (so they DO get used).
        # For BOTS mode the Mask will be 3D, to accommodate changes in masked
        # pixels per integration.
        log.info("Creating mask")
        Mask, nan_pix = create_mask(input_model, n_sigma)

        # Store the mask image in a model, if requested
        if save_mask:
            if len(Mask.shape) == 3:
                mask_model = datamodels.CubeModel(data=Mask)
            else:
                mask_model = datamodels.ImageModel(data=Mask)
        else:
            mask_model = None

    # If data are 3D, loop over integrations, applying correction to
    # each integration individually.
    log.info(f"Cleaning image {input_model.meta.filename}")
    if len(input_model.data.shape) == 3:

        # Check for 3D mask
        if len(Mask.shape) == 2:
            log.warning("Data are 3D, but Mask is 2D. Step will be skipped.")
            output_model.meta.cal_step.nsclean = 'SKIPPED'
            return output_model

        # Loop over integrations
        for i in range(Mask.shape[0]):
            log.debug(f" working on integration {i+1}")
            cleaner = NSClean(detector, Mask[i])
            image = np.float32(input_model.data[i])
            try:
                cleaned_image = cleaner.clean(image, buff=True)
                output_model.data[i] = cleaned_image
            except np.linalg.LinAlgError:
                log.warning("Error cleaning image; step will be skipped")
                output_model.meta.cal_step.nsclean = 'SKIPPED'
                return output_model

    # Clean a 2D image
    else:
        # Instantiate an NSClean object
        cleaner = NSClean(detector, Mask)
        image = np.float32(input_model.data)

        # Clean the image
        try:
            cleaned_image = cleaner.clean(image, buff=True)
            output_model.data = cleaned_image
        except np.linalg.LinAlgError:
            log.warning("Error cleaning image; step will be skipped")
            output_model.meta.cal_step.nsclean = 'SKIPPED'
            return output_model

    # Restore NaN's from original image
    output_model.data[nan_pix] = np.nan

    # Set completion status
    output_model.meta.cal_step.nsclean = 'COMPLETE'

    return output_model, mask_model
