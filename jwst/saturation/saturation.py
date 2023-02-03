
#  Module for 2d saturation
#
import logging
import numpy as np
from scipy.ndimage import binary_dilation

from stdatamodels.jwst.datamodels import dqflags

from ..lib import reffile_utils
from . import x_irs2

from stcal.saturation.saturation import flag_saturated_pixels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DONOTUSE = dqflags.pixel['DO_NOT_USE']
SATURATED = dqflags.pixel['SATURATED']
AD_FLOOR = dqflags.pixel['AD_FLOOR']
NO_SAT_CHECK = dqflags.pixel['NO_SAT_CHECK']
ATOD_LIMIT = 65535.  # Hard DN limit of 16-bit A-to-D converter


def flag_saturation(input_model, ref_model, n_pix_grow_sat):
    """
    Short Summary
    -------------
    Call function in stcal for flagging for saturated pixels.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    ref_model : `~jwst.datamodels.SaturationModel`
        Saturation reference file data model

    n_pix_grow_sat : int
        Number of layers of pixels adjacent to a saturated pixel to also flag
        as saturated (i.e '1' will flag the surrouding 8 pixels) to account for
        charge spilling.

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with saturation, A/D floor, and do not use flags set in
        the GROUPDQ array
    """

    data = input_model.data

    # Create the output model as a copy of the input
    output_model = input_model.copy()
    gdq = output_model.groupdq
    pdq = output_model.pixeldq

    zframe = input_model.zeroframe if input_model.meta.exposure.zero_frame else None

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

    gdq_new, pdq_new, zframe = flag_saturated_pixels(
        data, gdq, pdq, sat_thresh, sat_dq, ATOD_LIMIT, dqflags.pixel,
        n_pix_grow_sat=n_pix_grow_sat, zframe=zframe)

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = gdq_new

    # Save the NO_SAT_CHECK flags in the output PIXELDQ array
    output_model.pixeldq = pdq_new

    if zframe is not None:
        output_model.zeroframe = zframe

    return output_model


def irs2_flag_saturation(input_model, ref_model, n_pix_grow_sat):
    """
    Short Summary
    -------------
    For NIRSPEC IRS2 mode only, apply flagging for saturation based on threshold
    values stored in the saturation reference file and A/D floor based on
    testing for 0 DN values. For A/D floor flagged groups, the DO_NOT_USE flag
    is also set.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    ref_model : `~jwst.datamodels.SaturationModel`
        Saturation reference file data model

    n_pix_grow_sat : int
        Number of layers of pixels adjacent to a saturated pixel to also flag
        as saturated (i.e '1' will flag the surrouding 8 pixels) to account for
        charge spilling.

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with saturation, A/D floor, and do not use flags set in
        the GROUPDQ array
    """

    data = input_model.data
    nints = data.shape[0]
    ngroups = data.shape[1]
    detector = input_model.meta.instrument.detector

    # create a mask of the appropriate size
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
    # set the saturation check threshold to above the A-to-D converter limit,
    # so no pixels will ever be above that level and hence not get flagged.
    sat_thresh[np.bitwise_and(sat_dq, NO_SAT_CHECK) == NO_SAT_CHECK] = ATOD_LIMIT + 1

    # Also reset NaN values in the saturation threshold array to above
    # the A-to-D limit and flag them with NO_SAT_CHECK
    sat_dq[np.isnan(sat_thresh)] |= NO_SAT_CHECK
    sat_thresh[np.isnan(sat_thresh)] = ATOD_LIMIT + 1

    flagarray = np.zeros(data.shape[-2:], dtype=groupdq.dtype)
    flaglowarray = np.zeros(data.shape[-2:], dtype=groupdq.dtype)

    if input_model.meta.exposure.zero_frame:
        zflagarray = np.zeros(data.shape[-2:], dtype=groupdq.dtype)
        zflaglowarray = np.zeros(data.shape[-2:], dtype=groupdq.dtype)

    for ints in range(nints):
        for group in range(ngroups):
            # Update the 4D groupdq array with the saturation flag.
            sci_temp = x_irs2.from_irs2(data[ints, group, :, :],
                                        irs2_mask, detector)
            # check for saturation
            flag_temp = np.where(sci_temp >= sat_thresh, SATURATED, 0)
            # check for A/D floor
            flaglow_temp = np.where(sci_temp <= 0, AD_FLOOR | DONOTUSE, 0)

            # now, flag any pixels that border saturated pixels (not A/D floor pix)
            if n_pix_grow_sat > 0:
                flag_temp = adjacency_sat(flag_temp, SATURATED, n_pix_grow_sat)

            # Copy temps into flagarrays.
            x_irs2.to_irs2(flagarray, flag_temp, irs2_mask, detector)
            x_irs2.to_irs2(flaglowarray, flaglow_temp, irs2_mask, detector)

            # for saturation, the flag is set in the current plane
            # and all following planes.
            np.bitwise_or(groupdq[ints, group:, :, :], flagarray,
                          groupdq[ints, group:, :, :])
            # for A/D floor, the flag is only set of the current plane
            np.bitwise_or(groupdq[ints, group, :, :], flaglowarray,
                          groupdq[ints, group, :, :])

        # Process ZEROFRAME.  Instead of setting a ZEROFRAME DQ array, data
        # in the ZEROFRAME that is flagged will be set to 0.
        if input_model.meta.exposure.zero_frame:
            zplane = input_model.zeroframe[ints, :, :]
            zdq = np.zeros(groupdq.shape[-2:], dtype=groupdq.dtype)
            ztemp = x_irs2.from_irs2(zplane, irs2_mask, detector)

            zflag_temp = np.where(ztemp >= sat_thresh, SATURATED, 0)
            zflaglow_temp = np.where(ztemp <= 0, AD_FLOOR | DONOTUSE, 0)

            if n_pix_grow_sat > 0:
                zflag_temp = adjacency_sat(zflag_temp, SATURATED, n_pix_grow_sat)

            x_irs2.to_irs2(zflagarray, zflag_temp, irs2_mask, detector)
            x_irs2.to_irs2(zflaglowarray, zflaglow_temp, irs2_mask, detector)

            np.bitwise_or(zdq[:, :], zflagarray, zdq[:, :])
            np.bitwise_or(zdq[:, :], zflaglowarray, zdq[:, :])

            zplane[zdq != 0] = 0.
            output_model.zeroframe[ints, :, :] = zplane[:, :]
            del zdq

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = groupdq

    n_sat = np.any(np.any(np.bitwise_and(groupdq, SATURATED), axis=0), axis=0).sum()
    log.info(f'Detected {n_sat} saturated pixels')
    n_floor = np.any(np.any(np.bitwise_and(groupdq, AD_FLOOR), axis=0), axis=0).sum()
    log.info(f'Detected {n_floor} A/D floor pixels')

    # Save the NO_SAT_CHECK flags in the output PIXELDQ array
    pixeldq_temp = x_irs2.from_irs2(output_model.pixeldq, irs2_mask,
                                    detector)
    pixeldq_temp = np.bitwise_or(pixeldq_temp, sat_dq)
    x_irs2.to_irs2(output_model.pixeldq, pixeldq_temp, irs2_mask, detector)

    return output_model


def adjacency_sat(flag_temp, saturated, n_pix_grow_sat):
    """
    flag_temp : ndarray
        2D array of saturated groups.

    saturated : int
        Saturated flag.

    n_pix_grow_sat : int
        Number of layers of pixels adjacent to a saturated pixel to also flag
        as saturated (i.e '1' will flag the surrouding 8 pixels) to account for
        charge spilling.
    """
    only_sat = np.bitwise_and(flag_temp, saturated).astype(np.uint8)
    box_dim = (n_pix_grow_sat * 2) + 1
    struct = np.ones((box_dim, box_dim)).astype(bool)
    dialated = binary_dilation(only_sat, structure=struct).astype(only_sat.dtype)
    flag_temp = np.bitwise_or(flag_temp, (dialated * saturated))

    return flag_temp
