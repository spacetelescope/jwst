"""Pipeline utilities objects"""

import logging

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

from jwst.associations.lib.dms_base import TSO_EXP_TYPES


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def is_tso(model):
    """Is data Time Series Observation data?

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
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
    model : `~jwst.datamodels.JwstDataModel` or ndarray
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


def match_nans_and_flags(input_model):
    """Ensure data, error, variance, and DQ are marked consistently for invalid data.

    Invalid data is assumed to be any pixel set to NaN in any one of the
    data, error, or variance arrays, or else set to the DO_NOT_USE flag
    in the DQ array.

    The input model is updated in place with NaNs or DO_NOT_USE flags, as
    appropriate, at all invalid data locations.

    Parameters
    ----------
    input_model : DataModel
        Input model containing some combination of data, dq, err, var_rnoise,
        var_poisson, and var_flat extensions. These extensions must all have
        matching dimensions if present.
    """
    # Build up the invalid data flags from each available data extension.
    is_invalid = None
    data_shape = None
    nan_extensions = ['data', 'err', 'var_rnoise', 'var_poisson', 'var_flat']
    for extension in nan_extensions:
        if not hasattr(input_model, extension):
            continue
        data = getattr(input_model, extension)
        if is_invalid is None:
            is_invalid = np.isnan(data)
            data_shape = data.shape
        else:
            if data.shape != data_shape:
                log.warning(f"Mismatched data shapes; skipping invalid data "
                            f"updates for extension '{extension}'")
                continue
            is_invalid |= np.isnan(data)

    # Add in invalid flags from the DQ extension if present
    if hasattr(input_model, 'dq'):
        do_not_use = (input_model.dq & dqflags.pixel['DO_NOT_USE']).astype(bool)
        if is_invalid is None:
            is_invalid = do_not_use
            data_shape = input_model.dq.shape
        else:
            if input_model.dq.shape != data_shape:
                log.warning("Mismatched data shapes; skipping invalid data "
                            "updates for extension 'dq'")
            else:
                is_invalid |= do_not_use

    # Nothing to do if no extensions were found
    if is_invalid is None:
        return

    # Update all the data extensions
    for extension in nan_extensions:
        if not hasattr(input_model, extension):
            continue
        data = getattr(input_model, extension)
        if data.shape != data_shape:
            continue
        data[is_invalid] = np.nan

    # Update the DQ extension
    if input_model.dq.shape == data_shape:
        input_model.dq[is_invalid] |= dqflags.pixel['DO_NOT_USE']
