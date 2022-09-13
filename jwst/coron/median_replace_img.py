"""Replace bad pixels in the input images with the median of the surrounding pixels.
"""

import logging
import numpy as np
from jwst.datamodels import dqflags

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
        log.warning('Box for median filter is outside the data.')
        return 0.

    # Calculate the median value using only good pixels
    filtered_array = data_array[np.bitwise_and(dq_array, bad_bitvalue) == 0]
    median_value = np.median(filtered_array)

    if np.isnan(median_value):
        # If the median fails return 0
        log.debug('Median filter returned NaN; setting value to 0.')
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

        # Bad locations are defined as pixels set to NaN and/or pixels with
        # DQ flags contained in the list of bad dq bits
        bad_locations = np.where(np.bitwise_and(img_dq, bad_bitvalue) | np.isnan(img_int))

        # If it's MIRI coronagraphy, only use the bad locations that are in the
        # science aperture (i.e. not flagged as NON_SCIENCE). Set the others to 0
        # directly to avoid thousands of logging warnings.
        if img_model.meta.instrument.name == 'MIRI':
            bad_locations, non_science = separate_non_science_pixels(img_dq, bad_locations)
            # skip the median filter for non-science pixels
            img_int[non_science[0], non_science[1]] = 0

        # Fill the bad pixel values with the median of the data in the specified box region
        for i_pos in range(len(bad_locations[0])):
            # note: x and y are switched here but median_fill_value is
            # consistent with their usage here so it's all OK
            x_box_pos = bad_locations[0][i_pos]
            y_box_pos = bad_locations[1][i_pos]
            median_fill = median_fill_value(img_int, img_dq, box_size,
                                            bad_bitvalue, x_box_pos, y_box_pos)
            img_int[x_box_pos, y_box_pos] = median_fill

        img_model.data[nimage] = img_int

    return img_model


def separate_non_science_pixels(img_dq, bad_locations):
    """
    For the median filter, we don't care about the NON_SCIENCE pixels, but they
    produce a ton of warnings that clog up the alignment algorithm and make it
    run really slowly. In this function, we take all the bad pixels and pull out
    the ones with NON_SCIENCE flags so we can set them to 0 without running the
    median filter.

    Parameters
    ----------
    img_dq : ndarray
        input data quality array
    bad_locations : tuple
        2xN (row, col) tuple of flagged pixel indices

    Returns
    ------
    science_pixels : tuple
        2xN tuple of (row, col) flagged science pixels
    non_science_pixels : tuple
        2xN tuple of (row, col) flagged non_science pixels
    """

    def is_science(pix):
        # return True if pixel is for science, False if flagged NON_SCIENCE
        val = img_dq[pix[0], pix[1]]
        # check for the NON_SCIENCE flag
        flagged = np.bitwise_and(val, dqflags.pixel['NON_SCIENCE'])
        # if `flagged` is 0, it's a science pixel so return True.
        science_pixel = ~flagged.any()
        return science_pixel

    indexer = np.array(list(map(is_science, list(zip(*bad_locations)))))
    bad_locations = np.array(bad_locations)
    science_pixels = tuple(bad_locations[:, indexer])
    non_science_pixels = tuple(bad_locations[:, ~indexer])
    return science_pixels, non_science_pixels
