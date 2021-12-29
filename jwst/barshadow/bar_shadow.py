#
#  Module for calculating bar shadow correction for science data sets
#

import numpy as np
import logging
from gwcs import wcstools

from jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

SLITRATIO = 1.15     # Ratio of slit spacing to slit height


def do_correction(input_model, barshadow_model=None, inverse=False, source_type=None, correction_pars=None):
    """Do the Bar Shadow Correction

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        science data model to be corrected

    barshadow_model : `~jwst.datamodels.BarshadowModel`
        bar shadow data model from reference file

    inverse : boolean
        Invert the math operations used to apply the flat field.

    source_type : str or None
        Force processing using the specified source type.

    correction_pars : dict or None
        Correction parameters to use instead of recalculation.

    Returns
    -------
    output_model, corrections : `~jwst.datamodels.MultiSlitModel`, jwst.datamodels.DataModel
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
    log.debug('EXP_TYPE = %s' % exp_type)

    # Create output as a copy of the input science data model
    output_model = input_model.copy()

    # Loop over all the slits in the input model
    corrections = datamodels.MultiSlitModel()
    for slit_idx, slitlet in enumerate(output_model.slits):
        slitlet_number = slitlet.slitlet_id
        log.info('Working on slitlet %d' % slitlet_number)

        if correction_pars:
            correction = correction_pars.slits[slit_idx]
        else:
            correction = _calc_correction(slitlet, barshadow_model, source_type)
        corrections.slits.append(correction)

        # Apply the correction by dividing into the science and uncertainty arrays:
        #     var_poisson and var_rnoise are divided by correction**2,
        #     because they're variance, while err is standard deviation
        if not inverse:
            slitlet.data /= correction.data
        else:
            slitlet.data *= correction.data
        slitlet.err /= correction.data
        slitlet.var_poisson /= correction.data**2
        slitlet.var_rnoise /= correction.data**2
        if slitlet.var_flat is not None and np.size(slitlet.var_flat) > 0:
            slitlet.var_flat /= correction.data**2
        slitlet.barshadow = correction.data

    return output_model, corrections


def _calc_correction(slitlet, barshadow_model, source_type):
    """Calculate the barshadow correction for a slitlet

    Parameters
    ----------
    slitlet : jwst.datamodels.SlitModel
        The slitlet to calculate for.

    barshadow_model : `~jwst.datamodels.BarshadowModel`
        bar shadow data model from reference file

    source_type : str or None
        Force processing using the specified source type.

    Returns
    -------
    correction : jwst.datamodels.SlitModel
        The correction to be applied
    """
    slitlet_number = slitlet.slitlet_id

    # Create the pieces that are put together to make the barshadow model
    shutter_elements = create_shutter_elements(barshadow_model)
    w0 = barshadow_model.crval1
    wave_increment = barshadow_model.cdelt1
    y_increment = barshadow_model.cdelt2
    shutter_height = 1.0 / y_increment

    # The correction only applies to extended/uniform sources
    correction = datamodels.SlitModel(data=np.ones(slitlet.data.shape))
    if has_uniform_source(slitlet, source_type):
        shutter_status = slitlet.shutter_state
        if len(shutter_status) > 0:
            shadow = create_shadow(shutter_elements, shutter_status)

            # For each pixel in the slit subarray,
            # make a grid of indices for pixels in the subarray
            x, y = wcstools.grid_from_bounding_box(slitlet.meta.wcs.bounding_box, step=(1, 1))

            # Create the transformation from slit_frame to detector
            det2slit = slitlet.meta.wcs.get_transform('detector', 'slit_frame')

            # Use this transformation to calculate x, y, and wavelength
            xslit, yslit, wavelength = det2slit(x, y)

            # If the source position is off-center in the slit, renormalize the yslit
            # values so that it appears as if the source is centered, which is the appropriate
            # way to compute the shadow correction for extended/uniform sources (doesn't
            # depend on source location).
            if len(shutter_status) > 1:
                middle = (len(shutter_status) - 1) / 2.0
                src_loc = shutter_status.find('x')
                if src_loc != -1 and float(src_loc) != middle:
                    yslit -= np.nanmean(yslit)

            # The returned y values are scaled to where the slit height is 1
            # (i.e. a slit goes from -0.5 to 0.5).  The barshadow array is scaled
            # so that the separation between the slit centers is 1,
            # i.e. slit height + interslit bar
            yslit = yslit / SLITRATIO

            # Convert the Y and wavelength to a pixel location in the  bar shadow array;
            # the fiducial should always be at the center of the slitlet, regardless of
            # where the source is centered.
            index_of_fiducial = (len(shutter_status) - 1) / 2.0

            # The shutters go downwards, i.e. the first shutter in shutter_status corresponds to
            # the last in the shadow array.  So the center of the first shutter referred to in
            # shutter_status has an index of shadow.shape[0] - shutter_height.  Each subsequent
            # shutter center has an index shutter_height greater.
            index_of_fiducial_in_array = shadow.shape[0] - shutter_height * (1 + index_of_fiducial)
            yrow = index_of_fiducial_in_array + yslit * shutter_height
            wcol = (wavelength - w0) / wave_increment

            # Interpolate the bar shadow correction for non-Nan pixels
            correction.data = interpolate(yrow, wcol, shadow)
        else:
            log.info("Slitlet %d has zero length, correction skipped" % slitlet_number)
    else:
        log.info("Bar shadow correction skipped for slitlet %d (source not uniform)" % slitlet_number)

    return correction


def create_shutter_elements(barshadow_model):
    """Create the pieces that will be put together to make the barshadow
    array for the slitlets.  The pieces are:

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

       Parameters:

       barshadow_model: BarshadowModel object
           The barshadow model used to construct these pieces

       Returns:

       shutter_elements: dict
           Dictionary (specified above) with the pieces
    """
    shadow1x1 = barshadow_model.data1x1
    shadow1x3 = barshadow_model.data1x3
    shutter_elements = {}
    shutter_elements['first'] = create_first(shadow1x1)
    shutter_elements['open_open'] = create_open_open(shadow1x3)
    shutter_elements['open_closed'] = create_open_closed(shadow1x1)
    shutter_elements['closed_open'] = create_closed_open(shadow1x1)
    shutter_elements['closed_closed'] = create_closed_closed()
    shutter_elements['last'] = create_last(shadow1x1)

    return shutter_elements


def create_first(shadow1x1):
    """Create the first half shutter in the bar shadow array.
    Use rows 1-501 in the shadow1x1 array

    Parameters:

    shadow1x1: nddata array
        The 1x1 shadow array from the barshadow model

    Returns:

    The array to use as shutter_elements['first']
    """
    return shadow1x1[:501, :]


def create_open_open(shadow1x3):
    """Create the two half shutters obtained from two open shutters
    Use rows 1-501 in the shadow1x3 array

    Parameters:

    shadow1x3: nddata array
        The 1x3 shadow array from the barshadow model

    Returns:

    The array to use as shutter_elements['open_open']
    """
    return shadow1x3[:501, :]


def create_open_closed(shadow1x1):
    """Create the two half shutters obtained from one open and
    one closed shutter
    Use rows 501-1001 in the shadow1x1 array

    Parameters:

    shadow1x1: nddata array
        The 1x1 shadow array from the barshadow model

    Returns:

    The array to use as shutter_elements['open_closed']
    """
    return shadow1x1[500:, :]


def create_closed_open(shadow1x1):
    """Create the two half shutters obtained from one closed and
    one open shutter
    Use rows 1-501 in the shadow1x1 array

    Parameters:

    shadow1x1: nddata array
        The 1x1 shadow array from the barshadow model

    Returns:

    The array to use as shutter_elements['closed_open']
    """
    return shadow1x1[:501, :]


def create_closed_closed():
    """Create the two half shutters obtained from two closed shutters
    Uses 0.01 somewhat arbitrarily, although this case shouldn't occur
    very often, only when there are 2 or more consecutive closed shutters
    in a slitlet

    Parameters:

    None

    Returns:

    The array to use as shutter_elements['closed_closed']
    """
    return 0.01 * np.ones((501, 101))


def create_last(shadow1x1):
    """Create the last half shutter in the bar shadow array.
    Use rows 501-1001 in the shadow1x1 array

    Parameters:

    shadow1x1: nddata array
        The 1x1 shadow array from the barshadow model

    Returns:

    The array to use as shutter_elements['last']
    """
    return shadow1x1[501:, :]


def create_shadow(shutter_elements, shutter_status):
    """Create a bar shadow reference array on the fly from the shutter
    elements dictionary.

    Parameters:

    shutter_elements: dict
        The shutter elements dictionary

    shutter_status: string
        String describing the shutter status:
           0:  Closed
           1:  Open
           x:  Contains source

    Returns:

    shadow_array: nddata array
        The constructed bar shadow array

    """
    nshutters = len(shutter_status)
    shadow = create_empty_shadow_array(nshutters)
    first_row = 0
    shadow = add_first_half_shutter(shadow, shutter_elements['first'])
    first_row = first_row + shutter_elements['first'].shape[0] - 1
    last_shutter = 'open'
    for this_status in shutter_status[1:]:
        if this_status == '0':
            this_shutter = 'closed'
        else:
            this_shutter = 'open'
        shutter_pair = '_'.join((last_shutter, this_shutter))
        shadow = add_next_shutter(shadow, shutter_elements[shutter_pair], first_row)
        first_row = first_row + shutter_elements[shutter_pair].shape[0] - 1
        last_shutter = this_shutter
    shadow = add_last_half_shutter(shadow, shutter_elements['last'], first_row)
    return shadow


def create_empty_shadow_array(nshutters):
    """Create the empty bar shadow array.

    Parameters:

    nshutters: int
        The length of the slit in shutters

    Returns:

    empty_shadow: nddata_array
        The empty shadow array

    """
    #
    # Assume the reference files have a shape of 1001 rows by 101 columns
    # and go from -1 to +1 in Y
    nrows = nshutters * 500 + 500
    ncolumns = 101
    empty_shadow = np.zeros((nrows, ncolumns))
    return empty_shadow


def add_first_half_shutter(shadow, shadow_element):
    """Add the first half shutter to the shadow array

    Parameters:

    shadow: nddata array
        The bar shadow array.

    shadow_element: nddata array
        the shutter_elements['first'] array.  Should be 501 rows (Y) by 101 columns (wavelength)

    Returns:

    shadow: nddata array
        The bar shadow array with the first half shutter inserted

    """
    shadow[0:501, :] = shadow_element[:, :]
    return shadow


def add_next_shutter(shadow, shadow_element, first_row):
    """Add a single internal shutter and advance the last row by 501

    Parameters:

    shadow: nddata array
        The bar shadow array.

    shadow_element: nddata array
        the internal shutter element data array.

    first_row: int
        The first row to place the shutter element

    Returns:

    shadow: nddata array
        The bar shadow array with the double internal shutter inserted

    """
    #
    # Average the last row in the current bar shadow array with the first row of
    # the single internal shutter array
    shadow[first_row, :] = 0.5 * (shadow[first_row, :] + shadow_element[0, :])
    first_row = first_row + 1
    last_row = first_row + shadow_element.shape[0] - 1
    shadow[first_row:last_row, :] = shadow_element[1:, :]
    return shadow


def add_last_half_shutter(shadow, shadow_element, first_row):
    """Add the last half shutter from the 1x1 array.  The last half shutter
    is rows 501-1001.

    Parameters:

    shadow: nddata array
        The bar shadow array.

    shadow_element: nddata array
        the shadow_element array.

    first_row: int
        The first row to place the shadow element

    Returns:

    shadow: nddata array
        The bar shadow array with the last half shutter inserted

    """
    #
    # Average the last row in the current bar shadow array with the first row of
    # the shutter element array
    shadow[first_row, :] = 0.5 * (shadow[first_row, :] + shadow_element[0, :])
    first_row = first_row + 1
    last_row = first_row + shadow_element.shape[0] - 1
    shadow[first_row:last_row, :] = shadow_element[1:, :]
    return shadow


def interpolate(rows, columns, array, default=np.nan):
    """Interpolate row and column vectors in array

    Parameters:

    row: nddata array
        array of row indices

    column: nddata array
        array of column indices

    array: nddata array
        array to be interpolated

    default: number
        value to use in output array when input index is nan (default np.nan)

    Returns:

    correction: nddata array
        array of correction factors, or default when not calculated
    """
    nrows, ncolumns = rows.shape
    correction = np.ones((nrows, ncolumns))
    correction.fill(default)
    nrows_out, ncols_out = array.shape
    #
    # Extend the boundary of array by 1 row and column to handle end cases
    augmented_array = np.ones((nrows_out + 1, ncols_out + 1))
    augmented_array[:nrows_out, :ncols_out] = array
    augmented_array[nrows_out, :ncols_out] = array[nrows_out - 1, :]
    augmented_array[:nrows_out, ncols_out] = array[:, ncols_out - 1]
    augmented_array[nrows_out, ncols_out] = array[nrows_out - 1, ncols_out - 1]
    for row in range(nrows):
        for column in range(ncolumns):
            if ~np.isnan(rows[row, column]) and ~np.isnan(columns[row, column]):
                array_row = rows[row, column]
                array_column = columns[row, column]
                #
                # Deal with out-of-bounds pixels
                if array_row >= nrows_out:
                    array_row = nrows_out - 1
                if array_column >= ncols_out:
                    array_column = ncols_out - 1
                if array_row < 0:
                    array_row = 0
                if array_column < 0:
                    array_column = 0
                ix = int(array_column)
                iy = int(array_row)
                a11 = augmented_array[iy, ix]
                a12 = augmented_array[iy, ix + 1]
                a21 = augmented_array[iy + 1, ix]
                a22 = augmented_array[iy + 1, ix + 1]
                dx = array_column - ix
                dy = array_row - iy
                correction[row, column] = a11 * (1.0 - dx) * (1.0 - dy) + a12 * dx * (1.0 - dy) + \
                    a21 * (1.0 - dx) * dy + a22 * dx * dy
    return correction


def has_uniform_source(slitlet, force_type=None):
    """Determine whether the slitlet contains a uniform source

    Parameters:

    slitlet: slitlet object
        The slitlet being interrogated

    force_type : string or None
        Source type to force to and decide upon.

    Returns:

    answer: boolean
        True if the slitlet contains a uniform source
    """
    source_type = force_type if force_type else slitlet.source_type

    if source_type:
        # Assume extended, unless explicitly set to POINT
        if source_type.upper() == 'POINT':
            return False
        else:
            return True
    else:
        # If there's no source type info, default to EXTENDED
        log.info('SRCTYPE not set for slitlet %d; assuming EXTENDED' % slitlet.slitlet_id)
        return True
