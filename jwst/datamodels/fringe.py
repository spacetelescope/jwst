from .reference import ReferenceImageModel


__all__ = ['FringeModel']


class FringeModel(ReferenceImageModel):
    """
    A data model for 2D fringe correction images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    err : numpy array
        The error array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "fringe.schema.yaml"
