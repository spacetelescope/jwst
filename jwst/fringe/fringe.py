#
#  Module for applying fringe correction
#

import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, fringe_model):
    """
    Fringe-correct MIRI MRS data using a fringe model.

    Parameters
    ----------
    input_model : jwst.datamodels.IFUImageModel
        Input science data model to be fringe-corrected.

    fringe_model : jwst.datamodels.FringeModel
        Data model containing fringe.

    Returns
    -------
    output_model : jwst.datamodels.IFUImageModel
        Fringe-corrected science data model.
    """
    output_model = apply_fringe(input_model, fringe_model)
    output_model.meta.cal_step.fringe = "COMPLETE"

    return output_model


def apply_fringe(input_model, fringe):
    """
    Apply the fringe correction to the science data.

    The fringe correction corrects data and error arrays by dividing these
    arrays by the respective values in the reference array (for the reference
    pixels having good values), and updates the data quality array based on
    bad pixels in the reference array.

    Parameters
    ----------
    input_model : jwst.datamodels.IFUImageModel
        Input science data model to be fringe-corrected.

    fringe : jwst.datamodels.FringeModel
        Fringe reference file image model.

    Returns
    -------
    output_model : jwst.datamodel.IFUImageModel
        Input science data model which has been fringe-corrected.
    """
    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    # Fringe-correct data and error arrays, applying correction only
    # to pixels having good calibration values
    fringe_data = fringe.data
    fringe_data[np.isnan(fringe_data)] = 1.0
    inv_fringe_data_sq = 1 / (fringe_data * fringe_data)
    output_model.data /= fringe_data
    output_model.err /= fringe_data
    output_model.var_poisson *= inv_fringe_data_sq
    output_model.var_rnoise *= inv_fringe_data_sq
    if output_model.var_flat is not None and np.size(output_model.var_flat) > 0:
        output_model.var_flat *= inv_fringe_data_sq

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
    #   output_model.dq = np.bitwise_or(output_model.dq, dq_mask)
    log.info("DQ values in the reference file NOT used to update the output DQ.")

    return output_model
