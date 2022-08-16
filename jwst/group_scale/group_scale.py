#
#  Module for rescaling grouped data for NFRAME
#  not equal to a power of 2
#

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(model):
    """
    Short Summary
    -------------
    Rescales all groups in an exposure by FRMDIVSR/NFRAMES, to
    account for incorrect scaling when on-board frame averaging
    is done using a value of NFRAMES that is not a power of 2.
    The divisor for the on-board averaging is always a power of
    2, so if the number of frames coadded is not a power of 2,
    the result is incorrect and needs to be rescaled.

    Parameters
    ----------
    model: data model object
        science data to be corrected. Model is modified in place.
    """

    # Get the meta data values that we need
    nframes = model.meta.exposure.nframes
    frame_divisor = model.meta.exposure.frame_divisor
    if nframes is None or frame_divisor is None:
        log.warning('Necessary meta data not found')
        log.warning('Step will be skipped')
        model.meta.cal_step.group_scale = 'SKIPPED'
        return

    log.info('NFRAMES={}, FRMDIVSR={}'.format(nframes, frame_divisor))
    log.info('Rescaling all groups by {}/{}'.format(frame_divisor, nframes))

    # Apply the rescaling to the entire data array
    scale = float(frame_divisor) / nframes
    model.data *= scale
    model.meta.cal_step.group_scale = 'COMPLETE'

    return
