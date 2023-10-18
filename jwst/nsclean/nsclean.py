import logging
import numpy as np

from astropy.stats import sigma_clipped_stats

from jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.nsclean.lib import NSClean

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def mask_non_science(input_model, Mask):

    from .. assign_wcs import nirspec
    import gwcs

    log.info("Finding NON_SCIENCE pixels for an IFU image")

    # Initialize global DQ map
    dqmap = np.zeros_like(input_model.dq) + dqflags.pixel['NON_SCIENCE']

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
        # set non-NaN wavelength locations as good
        valid = ~np.isnan(wl)
        dq = dq[valid]
        x = x[valid]
        y = y[valid]
        dq[:] = 0

        # Copy dq for this slice into global dq map
        dqmap[y.astype(int), x.astype(int)] = dq

    # Now set all NON_SCIENCE locations in the mask to True
    nonsci_pix = np.where(dqmap & dqflags.pixel['NON_SCIENCE'])
    Mask[nonsci_pix] = True

    return Mask


def mask_slits(input_model, Mask):

    # Get the slit-to-msa frame transform from the WCS object
    slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')

    # Loop over the slits, marking all the pixels within each bounding
    # box as False (do not use) in the mask.
    # Note that for 3D masks, all planes will be set to the same value.
    for slit in slit2msa.slits:
        slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
        xlo, xhi, ylo, yhi = offset_wcs(slit_wcs)
        Mask[..., ylo:yhi, xlo:xhi] = False

    return Mask

    
def create_mask(input_model, n_sigma):
    """Create the pixel mask needed for setting which pixels to use
    for measuring 1/f noise
    """
    exptype = input_model.meta.exposure.type.lower()

    # Initialize mask
    if exptype == 'nrs_ifu':
        Mask = np.full(np.shape(input_model.dq), False)
    else:
        # Note: Mask will be 3D for BOTS mode data
        Mask = np.full(np.shape(input_model.dq), True)  # MOS and FS

    # If IFU, find NON_SCIENCE pixels and set to True in mask
    if exptype == 'nrs_ifu':
        Mask = mask_non_science(input_model, Mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if exptype in ['nrs_fixedslit', 'nrs_brightobj', 'nrs_msaspec']:
        Mask = mask_slits(input_model, Mask)

    # If IFU or MOS, set pixels affected by failed-open shutters to False
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        open_pix = np.where(input_model.dq & dqflags.pixel['MSA_FAILED_OPEN'])
        Mask[open_pix] = False

    # Reset NaN pixels and flag as don't use
    nan_pix = np.where(np.isnan(input_model.data))
    input_model.data[nan_pix] = 0
    Mask[nan_pix] = False

    # If IFU or MOS, don't use fixed-slit area pixels; uses hardwired indexes
    if exptype in ['nrs_ifu', 'nrs_msaspec']:
        Mask[922:1116, :] = False

    # Use left/right reference pixel columns (first and last 4). Can only be
    # applied to data that uses all 2048 columns of the detector.
    if Mask.shape[-1] == 2048:
        Mask[:, :5] = True
        Mask[:, -5:] = True  # keep one extra column on the right (always empty)

    # Flag outliers using sigma clipping stats.
    # For BOTS mode, which uses 3D data, loop over each integration separately.
    if len(input_model.data.shape) == 3:
        for i in range(input_model.data.shape[0]):
            _, median, sigma = sigma_clipped_stats(input_model.data[i], mask=~Mask[i], mask_value=0, sigma=5.0)
            outliers = np.where(input_model.data[i] > (median + n_sigma * sigma))
            Mask[i, outliers] = False
    else:
        _, median, sigma = sigma_clipped_stats(input_model.data, mask=~Mask, mask_value=0, sigma=5.0)
        outliers = np.where(input_model.data > (median + n_sigma * sigma))
        Mask[outliers] = False

    # Save the Mask, if desired
    if len(Mask.shape) == 3:
        mask_model = datamodels.CubeModel(data=Mask)
    else:
        mask_model = datamodels.ImageModel(data=Mask)
    #mask_model.save(input_model.meta.filename.replace('rate', 'mask'))
    mask_model.save(input_model.meta.filename.replace('msaflagopenstep', 'mask'))

    # Return the mask and the record of which pixels were NaN in the input;
    # it'll be needed later
    return Mask, nan_pix


def do_correction(input_model, n_sigma):

    """Execute the NSClean 1/f noise correction

    Parameters
    ----------
    input_model : data model object
        science data to be corrected

    n_sigma : float
        n-sigma rejection level for finding outliers

    Returns
    -------
    output_model : jwst.datamodel.JwstDataModel
        corrected data
    """

    detector = input_model.meta.instrument.detector.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')
    output_model = input_model.copy()

    # Create the pixel mask that'll be used to indicate which pixels
    # to include in the 1/f noise measurements. Basically, we're setting
    # all illuminated pixels to False, so that they do not get used, and
    # setting all unilluminated pixels to True (so they DO get used).
    # For BOTS mode the Mask will be 3D, to accommodate changes in masked
    # pixels per integration.
    log.info("Creating mask")
    Mask, nan_pix = create_mask(input_model, n_sigma)

    log.info(f"Cleaning image {input_model.meta.filename}")

    # If data are 3D, loop over integrations, applying correction to
    # each integration individually.
    if len(Mask.shape) == 3:
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

    return output_model
