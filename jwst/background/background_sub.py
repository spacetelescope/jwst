from __future__ import division

from .. import datamodels
from . import subtract_images

import numpy as np

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def background_sub(input_model, bkg_list):

    """
    Short Summary
    -------------
    Subtract the background signal from a JWST exposure by subtracting
    the average of one or more background exposures from it.

    Parameters
    ----------
    input_model: JWST data model
        input target exposure data model

    bkg_list: filename list
        list of background exposure file names

    Returns
    -------
    result: JWST data model
        background-subtracted target data model

    """

    # Compute the average of the background images associated with
    # the target exposure
    bkg_model = average_background(bkg_list)

    # Subtract the average background from the member
    log.debug(' subtracting avg bkg from %s', input_model.meta.filename)
    result = subtract_images.subtract(input_model, bkg_model)

    # Close the average background image and update the step status
    bkg_model.close()

    # We're done. Return the result.
    return result


def average_background(bkg_list):

    """
    Average multiple background exposures into a combined data model

    Parameters:
    -----------

    bkg_list: filename list
        List of background exposure file names

    Returns:
    --------

    avg_bkg: data model
        The averaged background exposure

    """

    avg_bkg = None

    # Loop over the images to be used as background
    for bkg_file in bkg_list:
        log.debug(' Accumulate bkg from %s', bkg_file)
        bkg_model = datamodels.ImageModel(bkg_file)

        # Initialize the avg_bkg model, if necessary
        if avg_bkg is None:
            avg_bkg = datamodels.ImageModel(bkg_model.data.shape)

        # Accumulate the data from this background image
        avg_bkg.data += bkg_model.data
        avg_bkg.err += bkg_model.err * bkg_model.err
        avg_bkg.dq = np.bitwise_or(avg_bkg.dq, bkg_model.dq)

        bkg_model.close()

    # Average the data in the accumulated background image
    num_bkg = len(bkg_list)
    avg_bkg.data = avg_bkg.data / num_bkg  # sci is normal average
    avg_bkg.err = np.sqrt(avg_bkg.err) / num_bkg  # err is uncertainty in the mean

    return avg_bkg
