"""Calculate bar shadow correction for science data sets."""

import numpy as np
import logging
from gwcs import wcstools
from scipy import ndimage
from stdatamodels.jwst import datamodels


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Fallback value for ratio of slit spacing to slit height
SLITRATIO = 1.15


def do_correction(
    input_model,
    barshadow_model=None,
    inverse=False,
    source_type=None,
    correction_pars=None,
):
    """
    Correct MSA data for bar shadows.

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        Science data model to be corrected.
    barshadow_model : `~jwst.datamodels.BarshadowModel`
        Bar shadow data model from reference file.
    inverse : bool
        Invert the math operations used to apply the flat field.
    source_type : str or None
        Force processing using the specified source type.
    correction_pars : dict or None
        Correction parameters to use instead of recalculation.

    Returns
    -------
    output_model, corrections : `~jwst.datamodels.MultiSlitModel`, jwst.datamodels.JwstDataModel
        Science data model with correction applied and barshadow extensions added,
        and a model of the correction arrays.
    """
    # Input is a MultiSlitModel science data model.
    # A MultislitModel has a member ".slits" that behaves like
    # a list of Slits, each of which has several data arrays and
    # keyword metadata items associated.
    #
    # At this point we're going to have to assume that a Slit is composed
    # of a set of slit[].nshutters in a vertical line.
    #
    # Reference file information is a 1x1 ref file and a 1x3 ref file.
    # Both the 1x1 and 1x3 ref files are 1001 pixels high and go from
    # -1 to 1 in their Y value WCS.

    exp_type = input_model.meta.exposure.type
    log.debug(f"EXP_TYPE = {exp_type}")

    # Create output as a copy of the input science data model
    output_model = input_model.copy()

    # Loop over all the slits in the input model
    corrections = datamodels.MultiSlitModel()
    for slit_idx, slitlet in enumerate(output_model.slits):
        slitlet_number = slitlet.slitlet_id
        log.info(f"Working on slitlet {slitlet_number}")

        if correction_pars:
            correction = correction_pars.slits[slit_idx]
        else:
            correction = _calc_correction(slitlet, barshadow_model, source_type)
        corrections.slits.append(correction)

        if correction is None:
            # For point sources, there is no real correction applied.
            # Record the correction status as False in this case.
            slitlet.barshadow_corrected = False

            # Store a blank barshadow image
            slitlet.barshadow = np.ones_like(slitlet.data)

            # No further processing needed
            continue

        # Otherwise, record the correction status as True.
        slitlet.barshadow_corrected = True

        # Apply the correction by dividing into the science and uncertainty arrays:
        #     var_poisson and var_rnoise are divided by correction**2,
        #     because they're variance, while err is standard deviation
        if not inverse:
            slitlet.data /= correction.data
            slitlet.err /= correction.data
            slitlet.var_poisson /= correction.data**2
            slitlet.var_rnoise /= correction.data**2
            if slitlet.var_flat is not None and np.size(slitlet.var_flat) > 0:
                slitlet.var_flat /= correction.data**2
        else:
            slitlet.data *= correction.data
            slitlet.err *= correction.data
            slitlet.var_poisson *= correction.data**2
            slitlet.var_rnoise *= correction.data**2
            if slitlet.var_flat is not None and np.size(slitlet.var_flat) > 0:
                slitlet.var_flat *= correction.data**2
        slitlet.barshadow = correction.data

    return output_model, corrections


def _calc_correction(slitlet, barshadow_model, source_type):
    """
    Calculate the barshadow correction for a slitlet.

    Parameters
    ----------
    slitlet : `~jwst.datamodels.SlitModel`
        The slitlet to calculate the correction for.

    barshadow_model : `~jwst.datamodels.BarshadowModel`
        Bar shadow data model from reference file.

    source_type : str or None
        Force processing using the specified source type.

    Returns
    -------
    correction : `~jwst.datamodels.SlitModel`
        The correction to be applied.
    """
    slitlet_number = slitlet.slitlet_id

    # Correction only applies to extended/uniform sources
    correction = None
    if not has_uniform_source(slitlet, source_type):
        log.info(f"Bar shadow correction skipped for slitlet {slitlet_number} (source not uniform)")
        return correction

    # No correction for zero length slitlets
    shutter_status = slitlet.shutter_state
    if len(shutter_status) == 0:
        log.info(f"Slitlet {slitlet_number} has zero length, correction skipped")
        return correction

    # Create the pieces that are put together to make the barshadow model
    shutter_elements = create_shutter_elements(barshadow_model)
    w0 = barshadow_model.crval1
    wave_increment = barshadow_model.cdelt1
    y_increment = barshadow_model.cdelt2
    shutter_height = 1.0 / y_increment

    shadow = create_shadow(shutter_elements, shutter_status)

    # For each pixel in the slit subarray,
    # make a grid of indices for pixels in the subarray
    x, y = wcstools.grid_from_bounding_box(slitlet.meta.wcs.bounding_box, step=(1, 1))

    # Create the transformation from slit_frame to detector
    det2slit = slitlet.meta.wcs.get_transform("detector", "slit_frame")

    # Use this transformation to calculate x, y, and wavelength
    xslit, yslit, wavelength = det2slit(x, y)

    # The returned y values are scaled to where the slit height is 1
    # (i.e. a slit goes from -0.5 to 0.5).  The barshadow array is scaled
    # so that the separation between the slit centers is 1,
    # i.e. slit height + interslit bar
    if slitlet.slit_yscale is None:
        log.warning(f"Slit height scale factor not found. Using default value {SLITRATIO}.")
        yslit = yslit / SLITRATIO
    else:
        yslit = yslit / slitlet.slit_yscale

    # Find the fiducial shutter to align the constructed shadow with the real array
    src_loc = shutter_status.find("x")
    index_of_fiducial = len(shutter_status) - src_loc

    # The shutters go downwards, i.e. the first shutter in shutter_status corresponds to
    # the last in the shadow array.  So the center of the first shutter referred to in
    # shutter_status has an index of shadow.shape[0] - shutter_height.  Each subsequent
    # shutter center has an index shutter_height greater.
    index_of_fiducial_in_array = shadow.shape[0] - shutter_height * index_of_fiducial
    yrow = index_of_fiducial_in_array + yslit * shutter_height
    wcol = (wavelength - w0) / wave_increment

    # Interpolate the bar shadow correction for non-Nan pixels
    correction = datamodels.SlitModel()
    correction.data = ndimage.map_coordinates(
        shadow, [yrow, wcol], cval=np.nan, order=1, mode="nearest"
    )

    return correction


def create_shutter_elements(barshadow_model):
    """
    Create half-shutter pieces for assembling a full barshadow array.

    The pieces are:

        1. shutter_elements['first']
           Goes from the bottom edge of the array (at 1 shutter width from
           the center of the first shutter) to the center of the first shutter
        2. shutter_elements['open_open']
            Goes from the center of an open shutter to the center of the next shutter,
            if that shutter is open
        3. shutter_elements['open_closed']
            Goes from the center of an open shutter to the center of the next shutter,
            if that shutter is closed
        4. shutter_elements['closed_open']
            Goes from the center of a closed shutter to the center of the next shutter,
            if that shutter is open
        5. shutter_elements['closed_closed']
            Goes from the center of a closed shutter to the center of the next shutter,
            if that shutter is also closed
        6. shutter_elements['last']
            Goes from the center of the last open shutter to the top edge of the shutter
            array (1 shutter width from the center of the last open shutter)

    Parameters
    ----------
    barshadow_model : BarshadowModel
        The barshadow model used to construct these pieces.

    Returns
    -------
    shutter_elements : dict
        Dictionary (specified above) with the pieces.
    """
    shadow1x1 = barshadow_model.data1x1
    shadow1x3 = barshadow_model.data1x3

    # Each shutter is assumed to be 1000 pixels tall, so
    # center elements are at row 500.
    # 1 pixel overlap is used at the boundaries, for averaging
    # corrections between adjoining sections
    # NOTE: the closed-closed element should not happen often, if ever:
    # it is rare to have two adjacent shutters closed within a slitlet
    shutter_elements = {
        "first": shadow1x1[:501, :],
        "open_open": shadow1x3[:501, :],
        "open_closed": shadow1x1[500:, :],
        "closed_open": shadow1x1[:501, :],
        "closed_closed": np.nanmin(shadow1x1) * np.ones((501, 101)),
        "last": shadow1x1[501:, :],
    }

    return shutter_elements


def create_shadow(shutter_elements, shutter_status):
    """
    Create a bar shadow reference array from the shutter elements.

    Parameters
    ----------
    shutter_elements : dict
        The shutter elements dictionary.
    shutter_status : str
        String describing the shutter status:
           0:  Closed
           1:  Open
           x:  Contains source

    Returns
    -------
    shadow_array : ndarray of float
        The constructed bar shadow array.
    """
    nshutters = len(shutter_status)
    shadow = create_empty_shadow_array(nshutters)
    first_row = 0
    shadow = add_first_half_shutter(shadow, shutter_elements["first"])
    first_row = first_row + shutter_elements["first"].shape[0] - 1
    last_shutter = "open"
    for this_status in shutter_status[1:]:
        if this_status == "0":
            this_shutter = "closed"
        else:
            this_shutter = "open"
        shutter_pair = "_".join((last_shutter, this_shutter))
        shadow = add_next_shutter(shadow, shutter_elements[shutter_pair], first_row)
        first_row = first_row + shutter_elements[shutter_pair].shape[0] - 1
        last_shutter = this_shutter
    shadow = add_last_half_shutter(shadow, shutter_elements["last"], first_row)
    return shadow


def create_empty_shadow_array(nshutters):
    """
    Create the empty bar shadow array.

    Parameters
    ----------
    nshutters : int
        The length of the slit in shutters.

    Returns
    -------
    empty_shadow : ndarray
        The empty shadow array.
    """
    # Assume the reference files have a shape of 1001 rows by 101 columns
    # and go from -1 to +1 in Y
    nrows = nshutters * 500 + 500
    ncolumns = 101
    empty_shadow = np.zeros((nrows, ncolumns))
    return empty_shadow


def add_first_half_shutter(shadow, shadow_element):
    """
    Add the first half shutter to the shadow array.

    Parameters
    ----------
    shadow : ndarray
        The bar shadow array.
    shadow_element : ndarray
        The shutter_elements['first'] array.  Should be 501 rows (Y)
        by 101 columns (wavelength).

    Returns
    -------
    shadow : ndarray
        The bar shadow array with the first half shutter inserted.
    """
    shadow[0:501, :] = shadow_element[:, :]
    return shadow


def add_next_shutter(shadow, shadow_element, first_row):
    """
    Add a single internal shutter and advance the last row by 501.

    Parameters
    ----------
    shadow : ndarray
        The bar shadow array.
    shadow_element : ndarray
        The internal shutter element data array.
    first_row : int
        The first row to place the shutter element.

    Returns
    -------
    shadow: ndarray
        The bar shadow array with the double internal shutter inserted.
    """
    # Average the last row in the current bar shadow array with the first row of
    # the single internal shutter array
    shadow[first_row, :] = 0.5 * (shadow[first_row, :] + shadow_element[0, :])
    first_row = first_row + 1
    last_row = first_row + shadow_element.shape[0] - 1
    shadow[first_row:last_row, :] = shadow_element[1:, :]
    return shadow


def add_last_half_shutter(shadow, shadow_element, first_row):
    """
    Add the last half shutter from the 1x1 array.

    The last half shutter is rows 501-1001.

    Parameters
    ----------
    shadow : nddata array
        The bar shadow array.
    shadow_element : nddata array
        The shadow_element array.
    first_row : int
        The first row to place the shadow element.

    Returns
    -------
    shadow: ndarray
        The bar shadow array with the last half shutter inserted.
    """
    #
    # Average the last row in the current bar shadow array with the first row of
    # the shutter element array
    shadow[first_row, :] = 0.5 * (shadow[first_row, :] + shadow_element[0, :])
    first_row = first_row + 1
    last_row = first_row + shadow_element.shape[0] - 1
    shadow[first_row:last_row, :] = shadow_element[1:, :]
    return shadow


def has_uniform_source(slitlet, force_type=None):
    """
    Determine whether the slitlet contains a uniform source.

    Parameters
    ----------
    slitlet : `~jwst.datamodels.SlitModel`
        The slitlet being interrogated.
    force_type : str or None
        Source type to force to and decide upon.

    Returns
    -------
    answer: bool
        True if the slitlet contains a uniform source.
    """
    source_type = force_type if force_type else slitlet.source_type

    if source_type:
        # Assume extended, unless explicitly set to POINT
        if source_type.upper() == "POINT":
            return False
        else:
            return True
    else:
        # If there's no source type info, default to EXTENDED
        log.info(f"SRCTYPE not set for slitlet {slitlet.slitlet_id}; assuming EXTENDED.")
        return True
