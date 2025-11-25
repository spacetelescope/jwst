import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.clean_flicker_noise import clean_flicker_noise as cfn

__all__ = ["correct_picture_frame"]

log = logging.getLogger(__name__)

LOWER_EDGE_REGION = (4, 75)
UPPER_EDGE_REGION = (1985, 2043)
CENTER_REGION = (76, 1984)


def _region_masks(background_mask):
    """
    Make region masks for center and edge data.

    Parameters
    ----------
    background_mask : ndarray
        Background mask image (1 = valid background pixel, 0 otherwise).

    Returns
    -------
    center_data : ndarray of bool
        Array matching background mask shape; True for background pixels
        in the center region.
    edge_data : ndarray of bool
        Array matching background mask shape; True for background pixels
        in the edge region.
    """
    y, x = np.mgrid[: background_mask.shape[0], : background_mask.shape[1]]

    # Background data in the center region
    center_data = (
        background_mask
        & (y >= CENTER_REGION[0])
        & (y < CENTER_REGION[1])
        & (x >= CENTER_REGION[0])
        & (x < CENTER_REGION[1])
    )

    # Background data in the edge regions, not including the reference pixels around the border
    out_of_range = (
        (y < LOWER_EDGE_REGION[0])
        | (y >= UPPER_EDGE_REGION[1])
        | (x < LOWER_EDGE_REGION[0])
        | (x >= UPPER_EDGE_REGION[1])
    )
    bottom = (y >= LOWER_EDGE_REGION[0]) & (y < LOWER_EDGE_REGION[1])
    top = (y >= UPPER_EDGE_REGION[0]) & (y < UPPER_EDGE_REGION[1])
    left = (x >= LOWER_EDGE_REGION[0]) & (x < LOWER_EDGE_REGION[1])
    right = (x >= UPPER_EDGE_REGION[0]) & (x < UPPER_EDGE_REGION[1])
    edge_data = background_mask & (bottom | top | left | right) & ~out_of_range

    return center_data, edge_data


def correct_picture_frame(
    input_model,
    pictureframe_model,
    mask_science_regions=True,
    n_sigma=2.0,
    input_dir="",
    save_mask=False,
    save_correction=False,
):
    """
    Correct full-frame NIRSpec data for the thermal picture frame effect.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel` \
                  or `~stdatamodels.jwst.datamodels.ImageModel` \
                  or `~stdatamodels.jwst.datamodels.CubeModel`
        Science model to be corrected. Updated in place.
    pictureframe_model : `~stdatamodels.jwst.datamodels.PictureFrameModel`
        Picture frame reference model.
    mask_science_regions : bool, optional
        Mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters.
    n_sigma : float, optional
        N-sigma rejection level for finding outliers.
    input_dir : str, optional
        Path to the input directory.  Used by sub-steps (e.g. assign_wcs
        for NIRSpec MOS data) to find auxiliary data.
    save_mask : bool, optional
        Switch to indicate whether the scene mask should be saved.
    save_correction : bool, optional
        Switch to indicate whether the scaled correction data should be saved.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel` \
                  or `~stdatamodels.jwst.datamodels.ImageModel` \
                  or `~stdatamodels.jwst.datamodels.CubeModel`
        Corrected data, updated in place.
    mask_model : `~stdatamodels.jwst.datamodels.ImageModel`
        Pixel mask to be saved or None.
    correction_model : `~stdatamodels.jwst.datamodels.RampModel` \
                       or `~stdatamodels.jwst.datamodels.ImageModel` \
                       or `~stdatamodels.jwst.datamodels.CubeModel`
        Correction model to be saved or None.
    status : {'COMPLETE', 'SKIPPED'}
        Completion status.  If errors were encountered, status = 'SKIPPED'
        and the output data is unchanged from the input data.  Otherwise,
        status = 'COMPLETE'.
    """
    # Set some default return values, to be updated later
    correction_model = None
    mask_model = None
    status = "SKIPPED"

    # Check for valid input: must be NIRSpec full frame
    instrument = str(input_model.meta.instrument.name).upper()
    subarray = str(input_model.meta.subarray.name).upper()
    if instrument != "NIRSPEC" or subarray != "FULL":
        log.warning("Picture frame correction is only applicable to NIRSpec full frame exposures.")
        log.warning("Processing will be skipped.")
        return input_model, mask_model, correction_model, status

    log.info("Applying picture frame correction")

    # Make a draft rate file if needed
    if isinstance(input_model, datamodels.RampModel):
        # Check for single group data: no processing is possible in this case
        if input_model.data.shape[1] == 1:
            log.warning("Picture frame correction cannot be performed for single-group data.")
            log.warning("Processing will be skipped.")
            return input_model, mask_model, correction_model, status

        image_model = cfn.make_rate(input_model)
        is_ramp = True
    else:
        # Input is already a rate or rateints file
        image_model = input_model
        is_ramp = False

    # Assign a WCS to the rate file and flag open MSA shutters
    try:
        image_model = cfn.post_process_rate(
            image_model,
            assign_wcs=mask_science_regions,
            msaflagopen=mask_science_regions,
            input_dir=input_dir,
        )
    except NoDataOnDetectorError:
        log.warning("WCS could not be assigned. Trying again with mask_science_regions=False.")
        mask_science_regions = False
        image_model = cfn.post_process_rate(
            image_model,
            assign_wcs=False,
            msaflagopen=False,
            input_dir=input_dir,
        )

    # Make a background mask from the draft rate file, blocking out science regions
    background_mask = cfn.create_mask(
        image_model, mask_science_regions=mask_science_regions, n_sigma=n_sigma, single_mask=True
    )
    if save_mask:
        mask_model = cfn.make_intermediate_model(image_model, background_mask)

    # Get the center and edge regions of the image, for computing statistics
    center_data, edge_data = _region_masks(background_mask)
    if not np.any(center_data) or not np.any(edge_data):
        log.error("No data to scale; skipping picture frame correction.")
        return input_model, mask_model, correction_model, status

    # Get the center and edge stats from the reference data
    slope = pictureframe_model.data
    median_center_ref = np.nanmedian(slope[center_data])
    median_edge_ref = np.nanmedian(slope[edge_data])
    log.debug(f"Median center reference: {median_center_ref}")
    log.debug(f"Median edge reference: {median_edge_ref}")

    # Keep a copy of the original input data to diff, if needed
    input_data_copy = None
    if save_correction:
        input_data_copy = input_model.data.copy()

    if is_ramp:
        # Skip processing the zero group, subtract it from each subsequent group
        zero_group = input_model.data[:, 0]
        image = input_model.data[:, 1:]
        image -= zero_group
    else:
        # Directly correct each image
        zero_group = 0.0
        image = input_model.data

    # Compute median values at the center and edges in each image
    median_center = np.nanmedian(image[..., center_data], axis=-1)
    median_edge = np.nanmedian(image[..., edge_data], axis=-1)
    log.debug(f"Median center: {median_center}")
    log.debug(f"Median edge: {median_edge}")

    # Scale reference data and subtract from science data
    scale = (median_center - median_edge) / (median_center_ref - median_edge_ref)
    difference = slope - median_edge_ref
    offset = median_edge

    # Scale and offset are both one value per group/integration; difference is a 2D image
    correction = scale[..., None, None] * difference[..., :, :] + offset[..., None, None]

    # Add back in the zero group and subtract the correction in place
    image += zero_group - correction

    # Save the correction data if desired
    if save_correction:
        correction_data = input_model.data - input_data_copy
        correction_model = cfn.make_intermediate_model(input_model, correction_data)

    # Correction successfully applied
    status = "COMPLETE"
    return input_model, mask_model, correction_model, status
