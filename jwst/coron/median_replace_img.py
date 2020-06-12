"""Replace bad pixels with the median of the surrounding pixel and median fill
   the input images.
 """
import logging
import numpy as np

from jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def median_fill_value(input_array, input_dq_array, bsize, xc, yc):
    """
        Arguments:
        ----------
        input_array : ndarray
            Input array to filter.
        input_dq_array : ndarray
            Input data quality array
        bsize : scalar
            box size of the data to extract
        xc: scalar
           x position of the data extraction
        xc: scalar
           y position of the data extraction
        """
    # set the half box size
    hbox = int(bsize/2)

    # Extract the region of interest for the data
    try:
        data_array = input_array[xc - hbox:xc + hbox, yc - hbox: yc + hbox]
        dq_array = input_dq_array[xc - hbox:xc + hbox, yc - hbox: yc + hbox]
    except IndexError:
        # If the box is outside the data return 0
        log.warning('Box for median filter is outside the data.')
        return 0.

    filtered_array = data_array[dq_array != dqflags.pixel['DO_NOT_USE']]
    median_value = np.median(filtered_array)

    if np.isnan(median_value):
        # If the median fails return 0
        log.warning('Median filter returned NaN setting value to 0.')
        median_value = 0.

    return median_value


def median_replace_img(img_model, box_size):
    """ Routine to replace any bad pixels with the median value of the surrounding
        pixels.
        Arguments:
        ----------
        input_array : image model
            Input array to filter.
        box_size : scalar
            box size for the median filter
    """

    n_ints, _, _ = img_model.data.shape
    for nimage in range(n_ints):
        img_int = img_model.data[nimage]
        img_dq = img_model.dq[nimage]
        # check to see if any of the pixels are flagged
        if np.count_nonzero(img_dq == dqflags.pixel['DO_NOT_USE']) > 0:
            bad_locations  = np.where(np.equal(img_dq, dqflags.pixel['DO_NOT_USE']))
            # fill the bad pixel values with the median of the data in a box region
            for i_pos in range(len(bad_locations[0])):
                x_box_pos = bad_locations[0][i_pos]
                y_box_pos = bad_locations[1][i_pos]
                median_fill = median_fill_value(img_int, img_dq, box_size, x_box_pos, y_box_pos)
                img_int[x_box_pos, y_box_pos] = median_fill

        img_model.data[nimage] = img_int

    return img_model
