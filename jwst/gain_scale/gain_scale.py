#
#  Module for rescaling data for non-standard gain value
#

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, gain_factor):
    """
    Short Summary
    -------------
    Rescales all integrations in an exposure by gain_factor, to
    account for non-standard detector gain settings. The SCI
    and ERR arrays are rescaled.

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    Returns
    -------
    output_model: data model object
        rescaled science data

    """

    # Create output as a copy of the input science data model
    output_model = input_model.copy()

    # Apply the rescaling to the entire data array
    log.info('Rescaling by {0}'.format(gain_factor))
    output_model.data *= gain_factor
    output_model.err *= gain_factor

    output_model.meta.exposure.gain_factor = gain_factor
    output_model.meta.cal_step.gain_scale = 'COMPLETE'

    return output_model
