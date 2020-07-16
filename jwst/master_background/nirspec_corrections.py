import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def correct_nrs_ifu_bkg(input):
    """Apply point source vs. uniform source pathloss adjustments
    to a NIRSpec IFU 2D master background array.

    Parameters
    ----------
    input : `~jwst.datamodels.IFUImageModel`
        The input background data.

    Returns
    -------
    input : `~jwst.datamodels.IFUIMAGEModel`
        An updated (in place) version of the input with the data
        replaced by the corrected 2D background.
    """

    log.info('Applying point source pathloss updates to IFU background')

    # Try to load the appropriate pathloss correction arrays
    try:
        pl_point = input.getarray_noinit('pathloss_point')
    except AttributeError:
        log.warning('Pathloss_point array not found in input')
        log.warning('Skipping pathloss background updates')
        return input

    try:
        pl_uniform = input.getarray_noinit('pathloss_uniform')
    except AttributeError:
        log.warning('Pathloss_uniform array not found in input')
        log.warning('Skipping pathloss background updates')
        return input

    # Apply the corrections
    input.data *= (pl_point / pl_uniform)

    return input
