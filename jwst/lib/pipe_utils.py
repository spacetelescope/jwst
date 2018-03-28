"""Pipeline utilities objects"""

from ..associations.lib.dms_base import TSO_EXP_TYPES
from ..datamodels import CubeModel

# Model types that typically represent TSO's
TSO_MODEL_TYPES = (CubeModel,)


def is_tso(model):
    """Is data Time Series Observation data?

    Parameters
    ----------
    model: datamodels.DataModel
        Data to check

    Returns
    -------
    is_tso: bool
       `True` if the model represents TSO data
    """
    is_tso = False

    # Check on exposure types
    try:
        is_tso = model.meta.exposure.type.lower() in TSO_EXP_TYPES
    except AttributeError:
        pass

    # Check on JWST-specific TSOVISIT flag
    try:
        is_tso = is_tso or model.meta.observation.tsovisit
    except AttributeError:
        pass

    # Check on DataModel types
    is_tso = is_tso or isinstance(model, TSO_MODEL_TYPES)

    return is_tso
