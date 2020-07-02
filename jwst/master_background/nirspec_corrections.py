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
    input.data *= (input.pathloss_point / input.pathloss_uniform)

    return input
