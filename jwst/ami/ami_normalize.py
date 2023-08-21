#
#  Module for normalizing the LG results for a science target by
#  the LG results for a reference target
#

import logging
from . import oifits

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def normalize_LG(target_model, reference_model):
    """
    Short Summary
    -------------
    Normalizes the LG results for a science target by the
    LG results for a reference target

    Parameters
    ----------
    target_model: AmiLgModel data model
        The target data to be normalized

    reference_model: AmiLgModel data model
        The reference data

    Returns
    -------
    output_model: AmiLgModel data model
        Normalized fringe data for the target
    """

    # Apply the normalizations to the target data
    # Initialize the calibration (normalization) class
    norm_model = oifits.CalibOifits(target_model, reference_model)
    output_model = norm_model.calibrate()

    # Return the normalized target model
    return output_model
