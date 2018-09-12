from . import image
from . import model_base

from .image import ImageModel


__all__ = ['IFUImageModel']


class IFUImageModel(ImageModel):
    """
    A data model for 2D IFU images.

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

    relsens2d: numpy array
        The relative sensitivity 2D array.
    """
    schema_url = "ifuimage.schema.yaml"
