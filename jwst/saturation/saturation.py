#
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

HUGE_NUM = 100000.

def do_correction(input_model, ref_model):
    """
    Short Summary
    -------------
    Execute all tasks for saturation, including using a saturation reference
    file.

    Parameters
    ----------
    input_model: data model object
        The input science data to be corrected

    ref_model: data model object
        Saturation reference file mode object

    Returns
    -------
    output_model: data model object
        object having GROUPDQ array saturation flags set
    """

    ramparr = input_model.data
    # Was IRS2 readout used?
    is_irs2_format = pipe_utils.is_irs2(input_model)
    if is_irs2_format:
        irs2_mask = x_irs2.make_mask(input_model)

   # Create the output model as a copy of the input
    output_model = input_model.copy()
    groupdq = output_model.groupdq

    # Extract subarray from reference file, if necessary
    if reffile_utils.ref_matches_sci(input_model, ref_model):
        satmask = ref_model.data
        dqmask = ref_model.dq
    else:
        log.info('Extracting reference file subarray to match science data')
        ref_sub_model = reffile_utils.get_subarray_model(input_model, ref_model)
        satmask = ref_sub_model.data.copy()
        dqmask = ref_sub_model.dq.copy()
        ref_sub_model.close()

    # For pixels flagged in reference file as NO_SAT_CHECK, set the dq mask
    #   and saturation mask
    wh_sat = np.bitwise_and(dqmask, dqflags.pixel['NO_SAT_CHECK'])
    dqmask[wh_sat == dqflags.pixel['NO_SAT_CHECK']] = dqflags.pixel['NO_SAT_CHECK']
    satmask[wh_sat == dqflags.pixel['NO_SAT_CHECK']] = HUGE_NUM
    # Correct saturation values for NaNs in the ref file
    correct_for_NaN(satmask, dqmask)

    dq_flag = dqflags.group['SATURATED']

    nints = ramparr.shape[0]
    ngroups = ramparr.shape[1]

    detector = input_model.meta.instrument.detector
    flagarray = np.zeros(ramparr.shape[-2:], dtype=groupdq.dtype)
    for ints in range(nints):
        for plane in range(ngroups):
            # Update the 4D groupdq array with the saturation flag. The
            # flag is set in the current plane and all following planes.
            if is_irs2_format:
                sci_temp = x_irs2.from_irs2(ramparr[ints, plane, :, :],
                                            irs2_mask, detector)
                flag_temp = np.where(sci_temp >= satmask, dq_flag, 0)
                # Copy flag_temp into flagarray.
                x_irs2.to_irs2(flagarray, flag_temp, irs2_mask, detector)
            else:
                flagarray[:, :] = np.where(ramparr[ints, plane, :, :] >= satmask,
                                          dq_flag, 0)
            np.bitwise_or(groupdq[ints, plane:, :, :], flagarray,
                          groupdq[ints, plane:, :, :])

    output_model.groupdq = groupdq
    if is_irs2_format:
        pixeldq_temp = x_irs2.from_irs2(output_model.pixeldq, irs2_mask,
                                        detector)
        pixeldq_temp = np.bitwise_or(pixeldq_temp, dqmask)
        x_irs2.to_irs2(output_model.pixeldq, pixeldq_temp, irs2_mask, detector)
    else:
        output_model.pixeldq = np.bitwise_or(output_model.pixeldq, dqmask)

    return output_model


def correct_for_NaN(satmask, dqmask):
    """
    Short Summary
    -------------
    If there are NaNs in the saturation values in the reference file, reset
    them to a very high value such that the comparison never results in a
    positive (saturated) result for the associated pixels in the science data.
    Also reset the associated dqmask values to indicate that, effectively,
    no saturation check will be done for those pixels.

    Parameters
    ----------
    satmask: 2-d array
        Subarray of saturation thresholds, from the saturation reference
        file.  This may be modified in-place.

    dqmask: ndarray, same shape as `satmask`
        The DQ array from the saturation reference file, used to update
        the PIXELDQ array in the output.  This may be modified in-place.
    """
    # If there are NaNs as the saturation values, update those values
    #     to ensure there will not be saturation.
    wh_nan = np.isnan(satmask)

    if np.any(wh_nan):
        satmask[wh_nan] = HUGE_NUM
        dqmask[wh_nan] |= dqflags.pixel['NO_SAT_CHECK']

        log.info("Unflagged pixels having saturation values set to NaN were"
                 " detected in the ref file; for those affected pixels no"
                 " saturation check will be made.")
