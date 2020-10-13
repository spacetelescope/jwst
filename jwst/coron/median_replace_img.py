"""Replace bad pixels in the input images with the median of the surrounding pixels.
"""

import logging
import numpy as np

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
    hbox = int(bsize/2)

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
        log.warning('Median filter returned NaN setting value to 0.')
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

        # Fill the bad pixel values with the median of the data in a box region
        for i_pos in range(len(bad_locations[0])):
            x_box_pos = bad_locations[0][i_pos]
            y_box_pos = bad_locations[1][i_pos]
            median_fill = median_fill_value(img_int, img_dq, box_size,
                                            bad_bitvalue, x_box_pos, y_box_pos)
            img_int[x_box_pos, y_box_pos] = median_fill

        img_model.data[nimage] = img_int

    return img_model
