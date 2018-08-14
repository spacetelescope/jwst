#
#  Module for calculating bar shadow correction for science data sets
#

import numpy as np
import logging
from .. import datamodels
from jwst.assign_wcs import nirspec
from gwcs import wcstools

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

SLITRATIO = 1.15     # Ratio of slit spacing to slit height

def do_correction(input_model, barshadow_model):
    """Do the Bar Shadow Correction

    Parameters
    ----------
    input_model: MultiSlitModel datamodel object
        science data model to be corrected

    barshadow_model: BarshadowModel datamodel object
        bar shadow datamodel from reference file

    Returns
    -------
    output_model: MultiSlitModel datamodel object
        Science datamodel with bar shadow extensions added

    """
#
# Input is a MultiSlitModel science data model 
# A MultislitModel has a member .slits that behaves like
# a list of Slits, each of which has several data arrays and
# keyword metadata items associated 
#
# At this point we're going to have to assume that a Slit is composed
# of a set of slit[].nshutters in a vertical line
#
# Reference file information is a 1x1 ref file and a 1x3 ref file
# Both the 1x1 and 1x3 ref files are 1001 pixel high and go from
# -1 to 1 in their Y value WCS.  
# 
# 
    exp_type = input_model.meta.exposure.type
    log.info(exp_type)
    #
    # Create the pieces that are put together to make the barshadow model
    shutter_elements = create_shutter_elements(barshadow_model)
    w0 = barshadow_model.crval1
    wave_increment = barshadow_model.cdelt1
    # For each slitlet
    for slitlet in input_model.slits:
        slitlet_number = slitlet.slitlet_id
        log.info('Working on slitlet %d' % slitlet_number)
        if has_uniform_source(slitlet):
            #
            # As Y increases, the pixel row number decreases, so the shutter_state is
            # 'upside down'
            shutter_status = slitlet.shutter_state[::-1]
            if len(shutter_status) > 0:
                shadow = create_shadow(shutter_elements, shutter_status)
                #
                # For each pixel in the slit subarray
                #   Make a grid of indices for pixels in the subarray
                x, y = wcstools.grid_from_bounding_box(slitlet.meta.wcs.bounding_box, step=(1,1))
                #   Create the transformation from slit_frame to detector
                det2slit = slitlet.meta.wcs.get_transform('detector', 'slit_frame')
                #   Use this transformation to calculate x, y and wavelength
                xslit, yslit, wavelength = det2slit(x, y)
                #   The returned y values are scaled to where the slit height is 1
                # (i.e. a slit goes from -0.5 to 0.5).  The barshadow array is scaled
                # so that the separation between the slit centers is 1, i.e. slit height
                # + interslit bar
                yslit = yslit / SLITRATIO
                #   Convert the Y and wavelength to a pixel location
                #   in the  bar shadow array
                index_of_fiducial = shutter_status.find('x')
                index_of_fiducial_in_array = 501 + index_of_fiducial*500
                yrow = index_of_fiducial_in_array - yslit*500.0
                wcol = (wavelength - w0)/wave_increment
                #   Interpolate the bar shadow correction for non-Nan pixels
                correction = interpolate(yrow, wcol, shadow)
                # Add the correction array and variance to the datamodel
                slitlet.barshadow = correction
            else:
                log.info("Slitlet %d has zero length, correction skipped" % slitlet_number)
                #
                # Put an array of ones in a correction extension
                slitlet.barshadow = np.ones(slitlet.data.shape)
        else:
            log.info("Bar shadow correction skipped for slitlet %d (source not uniform)" % slitlet_number)
            #
            # Put an array of ones in a correction extension
            slitlet.barshadow = np.ones(slitlet.data.shape)
    return input_model
    
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
#
#
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
    very often, only when there are 2 or more consecutuve closed shutters
    in a slitlet

    Parameters:

    None

    Returns:

    The array to use as shutter_elements['closed_closed']
    """
    return 0.01*np.ones(500)

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
    nrows = nshutters*500 + 500
    ncolumns = 101
    empty_shadow = np.zeros((nrows, ncolumns))
    return empty_shadow

def add_first_half_shutter(shadow, shadow_element):
    """Add the first half shutter to the shadow array

    Parameters:

    shadow: nddata array
        The bar shadow array.

    shadow_element: nddata array
        the shutter_elements['first'] array.  Should be 501 rows (Y) by 101 columns (wavelenth)

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
        The bar shadow array with the lastt half shutter inserted

    """
    #
    # Average the last row in the current bar shadow array with the first row of
    # the shutter element array
    shadow[first_row, :] = 0.5 * (shadow[first_row, :] + shadow_element[0, :])
    first_row = first_row + 1
    last_row = first_row + shadow_element.shape[0] - 1
    shadow[first_row:last_row, :] = shadow_element[1:, :]
    return shadow

def interpolate(rows, columns, array):
    """Interpolate row and column vectors in array

    Parameters:

    row: nddata array
        array of row indices

    column: nddata array
        array of column indices

    array: nddata array
        array to be interpolated
    """
    nrows, ncolumns = rows.shape
    correction = np.ones((nrows, ncolumns))
    nrows_out, ncols_out = array.shape
    #
    # Extend the boundary of array by 1 row and column to handle end cases
    augmented_array = np.ones((nrows_out+1, ncols_out+1))
    augmented_array[:nrows_out, :ncols_out] = array
    augmented_array[nrows_out,:ncols_out] = array[nrows_out-1, :]
    augmented_array[:nrows_out, ncols_out] = array[:, ncols_out-1]
    augmented_array[nrows_out, ncols_out] = array[nrows_out-1, ncols_out-1]
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
                a12 = augmented_array[iy, ix+1]
                a21 = augmented_array[iy+1, ix+1]
                a22 = augmented_array[iy+1, ix+1]
                dx = array_column - ix
                dy = array_row - iy
                correction[row, column] = a11*(1.0-dx)*(1.0-dy) + a12*dx*(1.0-dy) + \
                    a21*(1.0-dx)*dy + a22*dx*dy
    return correction

def has_uniform_source(slitlet):
    """Determine whether the slitlet contains a uniform source

    Parameters:

    slitlet: slitlet object
        The slitlet being interrogated

    Returns:

    answer: boolean
        True if the slitlet contains a uniform source
    """

    if slitlet.stellarity > 0.75:
        return False
    else:
        return True
