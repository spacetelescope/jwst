"""Pipeline utilities objects"""

import numpy as np

from ..associations.lib.dms_base import TSO_EXP_TYPES


def is_tso(model):
    """Is data Time Series Observation data?

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        Data to check

    Returns
    -------
    is_tso : bool
       `True` if the model represents TSO data
    """
    is_tso = False

    # Check on JWST-specific TSOVISIT flag
    try:
        is_tso = model.meta.visit.tsovisit
    except AttributeError:
        pass

    # Check on exposure types
    try:
        is_tso = is_tso or model.meta.exposure.type.lower() in TSO_EXP_TYPES
    except AttributeError:
        pass

    # Check on number of integrations
    try:
        if model.meta.exposure.nints is not None and model.meta.exposure.nints < 2:
            is_tso = False
    except AttributeError:
        pass

    # We've checked everything.
    return is_tso


def is_irs2(model):
    """Check whether the data are in IRS2 format.

    This currently assumes that only full-frame, near-infrared data can be
    taken using the IRS2 readout pattern.

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel` or ndarray
        Data to check

    Returns
    -------
    bool
       `True` if the data are in IRS2 format
    """

    if isinstance(model, np.ndarray):
        shape = model.shape
    else:
        try:
            shape = model.data.shape
        except AttributeError:
            return False

    max_length = 2048

    irs2_axis_length = max(shape[-1], shape[-2])

    if irs2_axis_length > max_length:
        return True
    else:
        return False
