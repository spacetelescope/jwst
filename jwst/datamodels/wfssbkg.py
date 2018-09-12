from .reference import ReferenceImageModel
from .dynamicdq import dynamic_mask


__all__ = ['WfssBkgModel']


class WfssBkgModel(ReferenceImageModel):
    """
    A data model for 2D WFSS master background reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  2-D.

    dq : numpy array
        The data quality array.  2-D.

    err : numpy array
        The error array.  2-D.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "wfssbkg.schema.yaml"
