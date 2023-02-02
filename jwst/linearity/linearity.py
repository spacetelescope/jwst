from stdatamodels.jwst.datamodels import dqflags
from ..lib import reffile_utils

import numpy as np
import logging

from stcal.linearity.linearity import linearity_correction

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, lin_model):

    # Create the output model as a copy of the input
    output_model = input_model.copy()
    zframe = None
    if output_model.meta.exposure.zero_frame:
        zframe = output_model.zeroframe

    # Get dq arrays
    pdq = input_model.pixeldq
    gdq = input_model.groupdq
    # If the input data does not have an expanded DQ array, create one
    if len(input_model.groupdq) == 0:
        gdq = (input_model.data * 0).astype(np.uint32)

    # Obtain linearity coefficients and dq array from reference file
    if reffile_utils.ref_matches_sci(input_model, lin_model):
        lin_coeffs = lin_model.coeffs
        lin_dq = lin_model.dq
    else:
        sub_lin_model = reffile_utils.get_subarray_model(input_model, lin_model)
        lin_coeffs = sub_lin_model.coeffs.copy()
        lin_dq = sub_lin_model.dq.copy()
        sub_lin_model.close()

    # Call linearity correction function in stcal
    new_data, new_pdq, new_zframe = linearity_correction(
        output_model.data, gdq, pdq, lin_coeffs, lin_dq, dqflags.pixel,
        zframe=zframe)

    output_model.data = new_data
    output_model.pixeldq = new_pdq
    if zframe is not None:
        output_model.zeroframe = new_zframe

    return output_model
