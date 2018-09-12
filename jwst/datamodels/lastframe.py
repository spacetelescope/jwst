from .reference import ReferenceImageModel


__all__ = ['LastFrameModel']


class LastFrameModel(ReferenceImageModel):
    """
    A data model for Last frame correction reference files.

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
    schema_url = "lastframe.schema.yaml"
