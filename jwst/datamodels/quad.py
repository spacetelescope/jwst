from .image import ImageModel

__all__ = ['QuadModel']


class QuadModel(ImageModel):
    """
    A data model for 4D image arrays.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  4-D.

    dq : numpy array
        The data quality array.  4-D.

    err : numpy array
        The error array.  4-D
    """
    schema_url = "quad.schema.yaml"
