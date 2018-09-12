from .reference import ReferenceImageModel


__all__ = ['FlatModel']


class FlatModel(ReferenceImageModel):
    """
    A data model for 2D flat-field images.

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
    schema_url = "flat.schema.yaml"
