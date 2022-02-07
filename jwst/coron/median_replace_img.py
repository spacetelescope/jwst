"""Replace bad pixels in the input images with the median of the surrounding pixels.
"""

import logging
import numpy as np

from pysiaf import Siaf

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def median_fill_value(input_array, input_dq_array, bsize, bad_bitvalue, xc, yc):
    """Calculates the median value of good pixels in a cutout of the input array.

    Parameters
    ----------
    input_array : ndarray
        Input array to filter

    input_dq_array : ndarray
        Input data quality array

    bsize : scalar
        Box size of the data to extract

    bad_bitvalue : int
        The sum of all of the DQ bit values to consider bad. Setting to 0
        will treat all pixels as good.

    xc : scalar
        x position of the data extraction

    yc : scalar
        y position of the data extraction

    Returns
    -------
    median_value : float
        The calculated median value
    """

    # Set the half box size
    hbox = int(bsize / 2)

    # Extract the region of interest for the data
    try:
        data_array = input_array[xc - hbox:xc + hbox + 1, yc - hbox: yc + hbox + 1]
        dq_array = input_dq_array[xc - hbox:xc + hbox + 1, yc - hbox: yc + hbox + 1]
    except IndexError:
        # If the box is outside the data return 0
        # log.warning('Box for median filter is outside the data.')
        return 0.

    # Calculate the median value using only good pixels
    filtered_array = data_array[np.bitwise_and(dq_array, bad_bitvalue) == 0]
    median_value = np.nanmedian(filtered_array)

    if np.isnan(median_value):
        # If the median fails return 0
        # log.warning('Median filter returned NaN setting value to 0.')
        median_value = 0.

    return median_value


def median_replace_img(img_model, box_size, bad_bitvalue):
    """Routine to replace any bad pixels with the median value of the surrounding
    pixels.

    Parameters
    ----------
    img_model : image model
        The input image model

    box_size : scalar
        Box size for the median filter

    bad_bitvalue : int
        The sum of all of the DQ bit values to consider bad. Setting to 0
        will treat all pixels as good.

    Returns
    -------
    img_model : image model

        The updated image model with the bad pixels replaced
    """

    n_ints, _, _ = img_model.data.shape
    for nimage in range(n_ints):
        img_int = img_model.data[nimage]
        img_dq = img_model.dq[nimage]
        bad_locations = np.where(np.bitwise_and(img_dq, bad_bitvalue))
        # if it's MIRI, only use the bad locations that are in the science aperture
        # set the others to 0 with no logging message
        bad_locations, zero_locations = filter_bad_locations(img_model, bad_locations)
        img_model.data[nimage, zero_locations[0], zero_locations[1]] = 0
        # Fill the bad pixel values with the median of the data in a box region
        for i_pos in range(len(bad_locations[0])):
            x_box_pos = bad_locations[0][i_pos]
            y_box_pos = bad_locations[1][i_pos]
            median_fill = median_fill_value(img_int, img_dq, box_size,
                                            bad_bitvalue, x_box_pos, y_box_pos)
            img_int[x_box_pos, y_box_pos] = median_fill

        img_model.data[nimage] = img_int

    return img_model


def check_coord(x, y, subarray, coron):
    """
    Check if a pixel whose coordinates are given in the subarray science frame
    is located within the coronagraph science aperture

    Parameters
    ----------
    x, y : pixel coordinate arrays
    subarray: relevant Siaf aperture object (.e.g MIRIM_MASK1550)
    coron: relevant Siaf aperture object (.e.g MIRIM_CORON1550)

    Output
    ------
    True if pixel is in the coronagraph science aperture, False otherwise
    """
    # convert from subarray to detector coordinates
    xy_det = subarray.convert(x, y, 'sci', 'det')
    # convert from detector to coronagraph science coordinates, because
    # the science coordinates aren't rotated
    xy_coron = coron.convert(*xy_det, 'det', 'sci')
    # get the limits
    xy_max = (np.max(coron.corners('sci'), axis=1)+0.5).astype(int)
    # make sure it's in the bounds, with some fancy array indexing for speed
    above_llim = np.array(xy_coron) >= 0 
    below_ulim = np.array(xy_coron) <= xy_max[:, None]
    is_inside  = np.concatenate([above_llim, below_ulim], axis=0).all(axis=0)
    return list(is_inside)


def filter_bad_locations(img_model, bad_locations):
    """
    Only use the bad locations that are in the science aperture. The rest, set to 0
    MIRI's coronagraphy readout subarrays have large regions of contiguous
    DO_NOT_USE pixels that can cause significant delays in the median filter
    if they are not handled intelligently.
    What we're going to do is only apply the median filter if the flagged pixel
    is in the MIRI coronagraphy science aperture; all flagged pixels outside
    science aperture simply get zeroed.

    Parameters
    ----------
    model : the datamodel

    Output
    ------
    bad_locations : locations to median filter
    ignore : locations to set to 0

    """
    # make sure you're only calling this on MIRI images
    name = img_model.meta.instrument.name
    coronagraph = img_model.meta.instrument.coronagraph.split("_")[0]
    try:
        assert(name == "MIRI")
        assert(coronagraph == "4QPM" or coronagraph == "LYOT")
    except:
        ignore = (np.array([], dtype=int), np.array([], dtype=int))
        return bad_locations, ignore
    coron_id = img_model.meta.instrument.coronagraph.split("_")[1]
    subarray_aper = Siaf("MIRI")[f'MIRIM_MASK{coron_id}']
    coron_aper = Siaf("MIRI")[f'MIRIM_CORON{coron_id}']
    is_in = check_coord(bad_locations[0], bad_locations[1], subarray_aper, coron_aper)
    # is_in = [check_coord(x, y, subarray_aper, coron_aper) for x,y in zip(*bad_locations)]
    still_bad = tuple(b[np.array(is_in)] for b in bad_locations)
    ignore = tuple(b[~np.array(is_in)] for b in bad_locations)
    return still_bad, ignore
