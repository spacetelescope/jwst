#
#  Module for normalizing the LG results for a science target by
#  the LG results for a reference target
#

import logging

from jwst.ami import oifits

log = logging.getLogger(__name__)

__all__ = ["normalize_lg"]


def normalize_lg(target_model, reference_model):
    """
    Normalize the LG results for a science target by the LG results for a reference target.

    Parameters
    ----------
    target_model : `~stdatamodels.jwst.datamodels.AmiOIModel`
        The target data to be normalized
    reference_model : `~stdatamodels.jwst.datamodels.AmiOIModel`
        The reference data

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.AmiOIModel`
        Normalized interferometric observables for the target
    """
    # Initialize the calibration (normalization) class and apply the normalizations
    norm_model = oifits.CalibOifits(target_model, reference_model)
    output_model = norm_model.calibrate()

    # Return the normalized target model
    return output_model
