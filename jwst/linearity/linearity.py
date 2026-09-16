import logging

import numpy as np
from stcal.linearity.linearity import linearity_correction
from stdatamodels.jwst.datamodels import dqflags

from jwst.lib import reffile_utils

log = logging.getLogger(__name__)

__all__ = ["do_correction"]


def do_correction(output_model, lin_model):
    """
    Apply the linearity correction to the data.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
        Science data to be corrected.

    lin_model : `~stdatamodels.jwst.datamodels.LinearityModel`
        Linearity reference file model.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel`
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
        gdq = np.zeros(output_model.data.shape, dtype=np.uint32)

    # Obtain linearity coefficients and dq array from reference file
    if reffile_utils.ref_matches_sci(output_model, lin_model):
        lin_coeffs = lin_model.coeffs
        lin_dq = lin_model.dq
        inv_coeffs = lin_model.inv_coeffs
    else:
        sub_lin_model = reffile_utils.get_subarray_model(output_model, lin_model)
        lin_coeffs = sub_lin_model.coeffs.copy()
        lin_dq = sub_lin_model.dq.copy()
        inv_coeffs = (
            sub_lin_model.inv_coeffs.copy() if sub_lin_model.inv_coeffs is not None else None
        )
        sub_lin_model.close()

    read_times_unset = (
        output_model.meta.exposure.read_times is None
        or len(output_model.meta.exposure.read_times) == 0
    )

    read_pattern = None
    if (
        inv_coeffs is not None  # Need to have inverse coefficients to correct
        and output_model.meta.exposure.nframes > 1  # Inv. correction for group-averaged data
        and read_times_unset  # If read times present, pattern is no longer applicable
    ):
        ngroups = output_model.meta.exposure.ngroups
        nframes = output_model.meta.exposure.nframes
        read_pattern = [
            [x + 1 + groupstart * nframes for x in range(nframes)] for groupstart in range(ngroups)
        ]
        log.info(
            "Providing inverse linearity coefficients and read_pattern to detailed "
            "linearity correction."
        )
    else:
        inv_coeffs = None  # Unset the inv_coeffs so stcal algo is not given this or read pattern

    # Call linearity correction function in stcal
    new_data, new_pdq, new_zframe = linearity_correction(
        output_model.data,
        gdq,
        pdq,
        lin_coeffs,
        lin_dq,
        dqflags.pixel,
        zframe=zframe,
        ilin_coeffs=inv_coeffs,
        read_pattern=read_pattern,
    )

    output_model.data = new_data
    output_model.pixeldq = new_pdq
    if zframe is not None:
        output_model.zeroframe = new_zframe

    return output_model
