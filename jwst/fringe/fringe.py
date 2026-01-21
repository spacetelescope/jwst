#
#  Module for applying fringe correction
#

import logging

import numpy as np

log = logging.getLogger(__name__)

__all__ = ["apply_fringe"]


def apply_fringe(input_model, fringe):
    """
    Apply the fringe correction to the science data.

    The fringe correction corrects data and error arrays by dividing these
    arrays by the respective values in the reference array (for the reference
    pixels having good values).

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        Input science data model to be fringe-corrected. Updated in place.

    fringe : `~stdatamodels.jwst.datamodels.FringeModel`
        Fringe reference file image model.

    Returns
    -------
    input_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
        Input science data model which has been fringe-corrected.
    """
    # Fringe-correct data and error arrays, applying correction only
    # to pixels having good calibration values
    fringe_data = fringe.data
    fringe_data[np.isnan(fringe_data)] = 1.0
    inv_fringe_data_sq = 1 / (fringe_data * fringe_data)
    input_model.data /= fringe_data
    input_model.err /= fringe_data
    input_model.var_poisson *= inv_fringe_data_sq
    input_model.var_rnoise *= inv_fringe_data_sq
    if input_model.var_flat is not None and np.size(input_model.var_flat) > 0:
        input_model.var_flat *= inv_fringe_data_sq

    # 05/22/14: For now, commenting out the following updating of the output
    # DQ based on the DQ of the reference file. This is done now because the
    # current DQ values in the ref file do not correspond to 'bad' data
    # values in the SCI array of the ref file.  Accordingly, this information
    # will be logged for now. This behavior may be changed later.
    ###
    # set DQ flag for bad pixels in the fringe
    #   dq_mask = fringe.dq * 0
    #   dq_mask[np.where(fringe.dq != 0)] = dqflags.pixel['DEAD']

    # update DQ array based on fringe's DQ values
    #   input_model.dq = np.bitwise_or(input_model.dq, dq_mask)
    log.debug("DQ values in the reference file NOT used to update the output DQ.")

    return input_model
