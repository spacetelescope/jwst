import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, gain_factor):
    """
    Short Summary
    -------------
    Rescales all integrations in an exposure by gain_factor, to
    account for non-standard detector gain settings. The SCI,
    ERR, and variance arrays are rescaled.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        Input datamodel to be corrected

    Returns
    -------
    output_model : `~jwst.datamodels.DataModel`
        Output datamodel with rescaled data

    """

    # Create output as a copy of the input science data model
    output_model = input_model.copy()

    # Apply the gain factor to the SCI and ERR arrays
    log.info('Rescaling by {0}'.format(gain_factor))
    output_model.data *= gain_factor
    output_model.err *= gain_factor

    # Apply the square of the gain factor to the variance arrays,
    # if they exist
    if (output_model.var_poisson is not None
            and np.size(output_model.var_poisson)) > 0:
        output_model.var_poisson *= gain_factor**2

    if (output_model.var_rnoise is not None
            and np.size(output_model.var_rnoise)) > 0:
        output_model.var_rnoise *= gain_factor**2

    # Set step status info
    output_model.meta.exposure.gain_factor = gain_factor
    output_model.meta.cal_step.gain_scale = 'COMPLETE'

    return output_model
