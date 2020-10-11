
#  Module for 2d saturation
#
import logging

from ..datamodels import dqflags
from ..lib import reffile_utils
from ..lib import pipe_utils
from . import x_irs2

import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

SATURATED = dqflags.pixel['SATURATED']
NO_SAT_CHECK = dqflags.pixel['NO_SAT_CHECK']
ATOD_LIMIT = 65535.  # Hard DN limit of 16-bit A-to-D converter


def do_correction(input_model, ref_model):
    """
    Short Summary
    -------------
    Apply saturation flagging based on threshold values stored in the
    saturation reference file.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    ref_model : `~jwst.datamodels.SaturationModel`
        Saturation reference file data model

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with saturation flags set in GROUPDQ array
    """

    ramp_array = input_model.data
    nints = ramp_array.shape[0]
    ngroups = ramp_array.shape[1]
    detector = input_model.meta.instrument.detector

    # If NIRSpec IRS2 readout mode was used, create a mask
    # of the appropriate size
    is_irs2_format = pipe_utils.is_irs2(input_model)
    if is_irs2_format:
        irs2_mask = x_irs2.make_mask(input_model)

    # Create the output model as a copy of the input
    output_model = input_model.copy()
    groupdq = output_model.groupdq

    # Extract subarray from saturation reference file, if necessary
    if reffile_utils.ref_matches_sci(input_model, ref_model):
        sat_thresh = ref_model.data
        sat_dq = ref_model.dq
    else:
        log.info('Extracting reference file subarray to match science data')
        ref_sub_model = reffile_utils.get_subarray_model(input_model, ref_model)
        sat_thresh = ref_sub_model.data.copy()
        sat_dq = ref_sub_model.dq.copy()
        ref_sub_model.close()

    # For pixels flagged in reference file as NO_SAT_CHECK,
    # set the saturation check threshold to the A-to-D converter limit.
    sat_thresh[np.bitwise_and(sat_dq, NO_SAT_CHECK) == NO_SAT_CHECK] = ATOD_LIMIT

    # Also reset NaN values in the saturation threshold array to the
    # A-to-D limit and flag them with NO_SAT_CHECK
    sat_dq[np.isnan(sat_thresh)] |= NO_SAT_CHECK
    sat_thresh[np.isnan(sat_thresh)] = ATOD_LIMIT

    # Loop over integrations and groups, checking for pixel values
    # that are above the saturation threshold
    flagarray = np.zeros(ramp_array.shape[-2:], dtype=groupdq.dtype)
    for ints in range(nints):
        for group in range(ngroups):
            # Update the 4D groupdq array with the saturation flag. The
            # flag is set in the current group and all following groups.
            if is_irs2_format:
                sci_temp = x_irs2.from_irs2(ramp_array[ints, group, :, :],
                                            irs2_mask,
                                            detector)
                flag_temp = np.where(sci_temp >= sat_thresh, SATURATED, 0)
                # Copy flag_temp into flagarray.
                x_irs2.to_irs2(flagarray, flag_temp, irs2_mask, detector)

            else:
                flagarray[:, :] = np.where(ramp_array[ints, group, :, :] >= sat_thresh,
                                           SATURATED, 0)
            np.bitwise_or(groupdq[ints, group:, :, :], flagarray,
                          groupdq[ints, group:, :, :])

    # Save the saturation flags in the output GROUPDQ array
    output_model.groupdq = groupdq

    n_sat = np.any(np.any(np.bitwise_and(groupdq, SATURATED), axis=0), axis=0).sum()
    log.info(f'Detected {n_sat} saturated pixels')

    # Save the NO_SAT_CHECK flags in the output PIXELDQ array
    if is_irs2_format:
        pixeldq_temp = x_irs2.from_irs2(output_model.pixeldq, irs2_mask,
                                        detector)
        pixeldq_temp = np.bitwise_or(pixeldq_temp, sat_dq)
        x_irs2.to_irs2(output_model.pixeldq, pixeldq_temp, irs2_mask, detector)
    else:
        output_model.pixeldq = np.bitwise_or(output_model.pixeldq, sat_dq)

    return output_model
