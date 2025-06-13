from stdatamodels.jwst.datamodels import dqflags
from jwst.lib import reffile_utils

import numpy as np
import logging

from stcal.linearity.linearity import linearity_correction

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(output_model, lin_model):
    """
    Apply the linearity correction to the data.

    Parameters
    ----------
    output_model : `~jwst.datamodels.RampModel`
        Science data to be corrected.

    lin_model : `~jwst.datamodels.LinearityModel`
        Linearity reference file model.

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Linearity corrected science data.
    """
    # Create the output model as a copy of the input
    zframe = None
    if output_model.meta.exposure.zero_frame:
        zframe = output_model.zeroframe

    # Get dq arrays
    pdq = output_model.pixeldq
    gdq = output_model.groupdq
    # If the input data does not have an expanded DQ array, create one
    if len(output_model.groupdq) == 0:
        gdq = (output_model.data * 0).astype(np.uint32)

    # Obtain linearity coefficients and dq array from reference file
    if reffile_utils.ref_matches_sci(output_model, lin_model):
        lin_coeffs = lin_model.coeffs
        lin_dq = lin_model.dq
    else:
        sub_lin_model = reffile_utils.get_subarray_model(output_model, lin_model)
        lin_coeffs = sub_lin_model.coeffs.copy()
        lin_dq = sub_lin_model.dq.copy()
        sub_lin_model.close()

    # Call linearity correction function in stcal
    new_data, new_pdq, new_zframe = linearity_correction(
        output_model.data, gdq, pdq, lin_coeffs, lin_dq, dqflags.pixel, zframe=zframe
    )

    output_model.data = new_data
    output_model.pixeldq = new_pdq
    if zframe is not None:
        output_model.zeroframe = new_zframe

    return output_model
