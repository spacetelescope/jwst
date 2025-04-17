#
#  Module for normalizing the LG results for a science target by
#  the LG results for a reference target
#

import logging
from . import oifits

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def normalize_lg(target_model, reference_model):
    """
    Normalize the LG results for a science target by the LG results for a reference target.

    Parameters
    ----------
    target_model : AmiOIModel data model
        The target data to be normalized
    reference_model : AmiOIModel data model
        The reference data

    Returns
    -------
    output_model : AmiOIModel data model
        Normalized interferometric observables for the target
    """
    # Initialize the calibration (normalization) class and apply the normalizations
    norm_model = oifits.CalibOifits(target_model, reference_model)
    output_model = norm_model.calibrate()

    # Return the normalized target model
    return output_model
